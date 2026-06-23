#include "pressure_poisson/pressure_poisson.hpp"

#include "core/particle_types.hpp"
#include "lsmps/lsmps_basis.hpp"
#include "lsmps/weight_function.hpp"

#include <petscksp.h>

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace lsmps {
namespace {

class PetscSession {
public:
    PetscSession() {
        PetscBool initialized = PETSC_FALSE;
        PetscInitialized(&initialized);
        if (!initialized) {
            int argc = 0;
            char** argv = nullptr;
            PetscInitialize(&argc, &argv, nullptr, nullptr);
            owns_session_ = true;
        }
    }

    ~PetscSession() {
        if (owns_session_) {
            PetscFinalize();
        }
    }

    PetscSession(const PetscSession&) = delete;
    PetscSession& operator=(const PetscSession&) = delete;

private:
    bool owns_session_ = false;
};

PetscSession& petscSession() {
    static PetscSession session;
    return session;
}

void checkPetsc(PetscErrorCode code, const char* operation) {
    if (code != PETSC_SUCCESS) {
        throw std::runtime_error(std::string("PETSc operation failed: ") + operation);
    }
}

std::vector<int> buildPressureDofMap(const ParticleSet& particles) {
    std::vector<int> dof(particles.size(), -1);
    int next = 0;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (particles.isFluid(i)) {
            dof[i] = next++;
        }
    }
    return dof;
}

int pressureDofCount(const std::vector<int>& dof) {
    int count = 0;
    for (const int value : dof) {
        if (value >= 0) {
            ++count;
        }
    }
    return count;
}

void accumulateVelocityRhs(
    const ParticleSet& particles,
    const ProvisionalVelocityResult& provisional,
    std::size_t i,
    std::size_t j,
    double support_radius,
    LsmpsKernelType kernel_type,
    LsmpsBasisVector& rhs_x,
    LsmpsBasisVector& rhs_y,
    LsmpsBasisVector& rhs_z) {
    const Vector3 offset = particles.positions()[j] - particles.positions()[i];
    const double weight = evaluateWeight(norm(offset), support_radius, kernel_type);
    const LsmpsBasisVector weighted_basis = weight * evaluateTypeABasis(offset, support_radius);
    const Vector3 velocity_delta = provisional.provisional_velocity[j] - provisional.provisional_velocity[i];
    rhs_x.noalias() += weighted_basis * velocity_delta.x;
    rhs_y.noalias() += weighted_basis * velocity_delta.y;
    rhs_z.noalias() += weighted_basis * velocity_delta.z;
}

double computeDivergence(
    const ParticleSet& particles,
    const TypedNeighborList& neighbors,
    const LsmpsInverseMatrix& regular,
    const ProvisionalVelocityResult& provisional,
    std::size_t i,
    double support_radius,
    LsmpsKernelType kernel_type) {
    LsmpsBasisVector rhs_x = LsmpsBasisVector::Zero();
    LsmpsBasisVector rhs_y = LsmpsBasisVector::Zero();
    LsmpsBasisVector rhs_z = LsmpsBasisVector::Zero();

    for (const std::size_t j : neighbors.fluid[i]) {
        accumulateVelocityRhs(particles, provisional, i, j, support_radius, kernel_type, rhs_x, rhs_y, rhs_z);
    }
    for (const std::size_t j : neighbors.wall[i]) {
        accumulateVelocityRhs(particles, provisional, i, j, support_radius, kernel_type, rhs_x, rhs_y, rhs_z);
    }

    const LsmpsBasisVector dx = regular.inverse_moment * rhs_x;
    const LsmpsBasisVector dy = regular.inverse_moment * rhs_y;
    const LsmpsBasisVector dz = regular.inverse_moment * rhs_z;
    return dx[0] / support_radius + dy[1] / support_radius + dz[2] / support_radius;
}

double laplacianCoefficient(
    const LsmpsInverseMatrix& pressure_matrix,
    const LsmpsBasisVector& basis,
    double weight,
    double support_radius) {
    const auto laplace_row =
        pressure_matrix.inverse_moment.row(3) + pressure_matrix.inverse_moment.row(4) + pressure_matrix.inverse_moment.row(5);
    return weight * (laplace_row * basis)(0) / (support_radius * support_radius);
}

double wallSourceCoefficient(
    const LsmpsInverseMatrix& pressure_matrix,
    const LsmpsBasisVector& q_basis,
    double weight,
    double support_radius) {
    const auto laplace_row =
        pressure_matrix.inverse_moment.row(3) + pressure_matrix.inverse_moment.row(4) + pressure_matrix.inverse_moment.row(5);
    return weight * (laplace_row * q_basis)(0) / support_radius;
}

bool isDirichletPressureParticle(const ParticleSet& particles, std::size_t i) {
    return particles.fluidStates()[i] == FluidParticleState::FreeSurface;
}

void setMatrixValue(Mat matrix, PetscInt row, PetscInt col, PetscScalar value) {
    checkPetsc(MatSetValues(matrix, 1, &row, 1, &col, &value, ADD_VALUES), "MatSetValues");
}

void setMatrixRow(
    Mat matrix,
    PetscInt row,
    const std::vector<PetscInt>& columns,
    const std::vector<PetscScalar>& values) {
    if (columns.empty()) {
        return;
    }
    checkPetsc(
        MatSetValues(
            matrix,
            1,
            &row,
            static_cast<PetscInt>(columns.size()),
            columns.data(),
            values.data(),
            ADD_VALUES),
        "MatSetValues row");
}

void setVectorValue(Vec vector, PetscInt row, PetscScalar value, InsertMode mode = INSERT_VALUES) {
    checkPetsc(VecSetValues(vector, 1, &row, &value, mode), "VecSetValues");
}

}  // namespace

PressurePoissonResult PressurePoissonAssembler::solve(
    const ParticleSet& particles,
    const TypedNeighborList& neighbors,
    const LsmpsMatrixSet& matrices,
    const ProvisionalVelocityResult& provisional,
    const SimulationConfig& config) const {
    if (!provisional.computed || provisional.provisional_velocity.size() != particles.size()) {
        throw std::runtime_error("PressurePoissonAssembler requires computed provisional velocities");
    }

    (void)petscSession();

    PressurePoissonResult result;
    result.pressure = particles.pressures();
    result.diagnostics.pressure_dof = buildPressureDofMap(particles);
    const int dof_count = pressureDofCount(result.diagnostics.pressure_dof);
    result.diagnostics.rhs.assign(particles.size(), 0.0);
    result.diagnostics.divergence.assign(particles.size(), 0.0);
    result.diagnostics.wall_neumann_source.assign(particles.size(), 0.0);

    Mat matrix = nullptr;
    Vec rhs = nullptr;
    Vec solution = nullptr;
    KSP ksp = nullptr;

    checkPetsc(MatCreate(PETSC_COMM_WORLD, &matrix), "MatCreate");
    checkPetsc(MatSetSizes(matrix, PETSC_DECIDE, PETSC_DECIDE, dof_count, dof_count), "MatSetSizes");
    checkPetsc(MatSetFromOptions(matrix), "MatSetFromOptions");
    checkPetsc(MatSeqAIJSetPreallocation(matrix, 96, nullptr), "MatSeqAIJSetPreallocation");
    checkPetsc(MatMPIAIJSetPreallocation(matrix, 96, nullptr, 96, nullptr), "MatMPIAIJSetPreallocation");
    checkPetsc(MatSetUp(matrix), "MatSetUp");
    checkPetsc(MatSetOption(matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE), "MatSetOption new nonzero");

    checkPetsc(VecCreate(PETSC_COMM_WORLD, &rhs), "VecCreate rhs");
    checkPetsc(VecSetSizes(rhs, PETSC_DECIDE, dof_count), "VecSetSizes rhs");
    checkPetsc(VecSetFromOptions(rhs), "VecSetFromOptions rhs");
    checkPetsc(VecDuplicate(rhs, &solution), "VecDuplicate solution");

    const double re = config.geometry.support_radius;
    std::vector<PetscInt> row_columns;
    std::vector<PetscScalar> row_values;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (!particles.isFluid(i)) {
            continue;
        }

        const PetscInt row = result.diagnostics.pressure_dof[i];
        if (isDirichletPressureParticle(particles, i)) {
            setMatrixValue(matrix, row, row, 1.0);
            setVectorValue(rhs, row, 0.0);
            result.diagnostics.rhs[i] = 0.0;
            continue;
        }

        const LsmpsInverseMatrix& regular = matrices.particles[i].regular;
        const LsmpsInverseMatrix& pressure_matrix = matrices.particles[i].pressure_neumann;
        if ((regular.status != LsmpsMatrixStatus::Valid && regular.status != LsmpsMatrixStatus::IllConditioned) ||
            (pressure_matrix.status != LsmpsMatrixStatus::Valid && pressure_matrix.status != LsmpsMatrixStatus::IllConditioned)) {
            setMatrixValue(matrix, row, row, 1.0);
            setVectorValue(rhs, row, result.pressure[i]);
            result.diagnostics.rhs[i] = result.pressure[i];
            continue;
        }

        const double divergence =
            computeDivergence(particles, neighbors, regular, provisional, i, re, config.lsmps.kernel_type);
        result.diagnostics.divergence[i] = divergence;

        double diagonal = 0.0;
        row_columns.clear();
        row_values.clear();
        row_columns.reserve(neighbors.fluid[i].size() + 1);
        row_values.reserve(neighbors.fluid[i].size() + 1);
        const Vector3 position_i = particles.positions()[i];
        for (const std::size_t j : neighbors.fluid[i]) {
            if (!particles.isFluid(j)) {
                continue;
            }
            const Vector3 offset = particles.positions()[j] - position_i;
            const double weight = evaluateWeight(norm(offset), re, config.lsmps.kernel_type);
            const double coefficient =
                laplacianCoefficient(pressure_matrix, evaluateTypeABasis(offset, re), weight, re);
            const PetscInt col = result.diagnostics.pressure_dof[j];
            if (col >= 0) {
                row_columns.push_back(col);
                row_values.push_back(coefficient);
                diagonal -= coefficient;
            }
        }

        double wall_source = 0.0;
        for (const std::size_t j : neighbors.wall[i]) {
            const Vector3 offset = particles.positions()[j] - position_i;
            const double weight = evaluateWeight(norm(offset), re, config.lsmps.kernel_type);
            const LsmpsBasisVector q_basis = evaluateTypeANeumannBasis(offset, particles.wallNormals()[j], re);
            const double neumann = config.physical.density * dot(config.physical.gravity, particles.wallNormals()[j]);
            wall_source += wallSourceCoefficient(pressure_matrix, q_basis, weight, re) * neumann;
        }
        result.diagnostics.wall_neumann_source[i] = wall_source;

        row_columns.push_back(row);
        row_values.push_back(diagonal);
        setMatrixRow(matrix, row, row_columns, row_values);
        const double rhs_value = config.physical.density / config.time.dt * divergence - wall_source;
        setVectorValue(rhs, row, rhs_value);
        result.diagnostics.rhs[i] = rhs_value;
    }

    checkPetsc(MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY), "MatAssemblyBegin");
    checkPetsc(MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY), "MatAssemblyEnd");
    checkPetsc(VecAssemblyBegin(rhs), "VecAssemblyBegin rhs");
    checkPetsc(VecAssemblyEnd(rhs), "VecAssemblyEnd rhs");

    checkPetsc(KSPCreate(PETSC_COMM_WORLD, &ksp), "KSPCreate");
    checkPetsc(KSPSetOperators(ksp, matrix, matrix), "KSPSetOperators");
    checkPetsc(KSPSetType(ksp, KSPGMRES), "KSPSetType");
    PC pc = nullptr;
    checkPetsc(KSPGetPC(ksp, &pc), "KSPGetPC");
    checkPetsc(PCSetType(pc, PCILU), "PCSetType");
    checkPetsc(KSPSetTolerances(ksp, config.linear_solver.tolerance, PETSC_DEFAULT, PETSC_DEFAULT, static_cast<PetscInt>(config.linear_solver.max_iterations)), "KSPSetTolerances");
    checkPetsc(KSPSetFromOptions(ksp), "KSPSetFromOptions");
    checkPetsc(KSPSolve(ksp, rhs, solution), "KSPSolve");

    PetscInt iterations = 0;
    PetscReal residual_norm = 0.0;
    KSPConvergedReason reason;
    checkPetsc(KSPGetIterationNumber(ksp, &iterations), "KSPGetIterationNumber");
    checkPetsc(KSPGetResidualNorm(ksp, &residual_norm), "KSPGetResidualNorm");
    checkPetsc(KSPGetConvergedReason(ksp, &reason), "KSPGetConvergedReason");
    result.solve_info.iterations = static_cast<int>(iterations);
    result.solve_info.final_residual_norm = static_cast<double>(residual_norm);
    result.solve_info.converged_reason = static_cast<int>(reason);
    result.solve_info.converged = reason > 0;
    result.solved = result.solve_info.converged;

    const PetscScalar* values = nullptr;
    checkPetsc(VecGetArrayRead(solution, &values), "VecGetArrayRead");
    for (std::size_t i = 0; i < particles.size(); ++i) {
        const int dof = result.diagnostics.pressure_dof[i];
        if (dof >= 0) {
            result.pressure[i] = PetscRealPart(values[dof]);
        }
    }
    checkPetsc(VecRestoreArrayRead(solution, &values), "VecRestoreArrayRead");

    checkPetsc(KSPDestroy(&ksp), "KSPDestroy");
    checkPetsc(VecDestroy(&solution), "VecDestroy solution");
    checkPetsc(VecDestroy(&rhs), "VecDestroy rhs");
    checkPetsc(MatDestroy(&matrix), "MatDestroy");

    return result;
}

}  // namespace lsmps
