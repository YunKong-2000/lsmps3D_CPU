# PPE 构造与求解模块设计

本文档确定压力泊松方程 PPE 的线性代数库、稀疏矩阵格式、迭代求解器选择和第一版实现路径。

## 1. 需求

PPE 线性系统来自 LSMPS type-A 压力拉普拉斯离散：

```text
A p = b
```

当前需求：

- 矩阵规模会随粒子数快速增大，必须使用稀疏矩阵。
- LSMPS PPE 系数矩阵通常不对称，不能默认使用 CG。
- 只能使用迭代方法。
- 尽可能快和准确。
- 支持多核 CPU。
- 后续需要可视化残差、迭代次数、收敛状态和每粒子压力。

## 2. 线性代数库选择

第一版选择 PETSc。

理由：

- PETSc 面向大型稀疏线性系统，核心 `KSP` 接口统一支持多种迭代求解器。
- 支持非对称系统常用的 GMRES、FGMRES、BiCGStab 等 Krylov 方法。
- 支持 AIJ/CSR 稀疏矩阵格式，适合按粒子邻居逐行装配 PPE。
- 支持 MPI 并行。当前 Homebrew 安装的 PETSc 3.24.2 可通过 OpenMPI 使用多进程并行。
- 求解器、预条件器、重启步数、容差、监控输出等可以通过 PETSc options 在运行时切换，便于后续调参。
- 后续若需要更强预条件器，可以切换到带 HYPRE 的 PETSc 构建并使用 BoomerAMG。

当前本机可用配置：

```text
PETSc version: 3.24.2
pkg-config: PETSc
compile/link flags:
    -I/opt/homebrew/opt/petsc/include -L/opt/homebrew/opt/petsc/lib -lpetsc
MPI: OpenMPI 可用
HYPRE: 当前 Homebrew PETSc 未启用
```

注意：`pkg-config PETSc` 当前只给出 PETSc 自身 include/lib 参数，不包含 MPI include 路径。直接使用普通 `c++` 编译会找不到 `mpi.h`。项目 CMake 接入 PETSc 时必须同时接入 MPI，或使用 MPI compiler wrapper。

不选择其他库作为第一版默认方案：

- Eigen：适合局部 9x9 矩阵和小规模稀疏验证，但大型非对称 PPE 的并行迭代求解与预条件能力不足。
- SuiteSparse：以直接法为主，不符合“只能使用迭代方法”的要求。
- OpenBLAS：主要是 BLAS/LAPACK 后端，不提供完整稀疏 Krylov 求解框架。
- 自写 CSR + Krylov：维护成本高，预条件、并行、诊断和鲁棒性都不如 PETSc。

## 3. 求解器策略

第一版默认求解器建议：

```text
KSP: GMRES 或 FGMRES
PC:  ILU（单进程）/ BJACOBI + ILU（MPI 多进程）
matrix: AIJ
```

默认推荐：

```text
ksp_type = gmres
pc_type = ilu
restart = 50 或 100
relative_tolerance = config.linear_solver.tolerance
max_iterations = config.linear_solver.max_iterations
```

说明：

- GMRES 是非对称线性系统的稳健第一选择。
- FGMRES 可作为后续默认候选，适合可变预条件器或更复杂预条件策略。
- BiCGStab 内存占用更低，可作为大规模算例中的 fallback，但收敛稳定性通常不如 GMRES 可控。
- 当前 PETSc 未带 HYPRE，因此第一版不默认使用 AMG。

运行时应允许通过 PETSc options 覆盖默认值，例如：

```text
-ksp_type gmres
-ksp_gmres_restart 100
-pc_type bjacobi
-sub_pc_type ilu
-ksp_rtol 1e-10
-ksp_max_it 1000
-ksp_monitor
```

## 4. PPE 矩阵格式

PETSc 侧使用 AIJ 稀疏矩阵：

```text
Mat A
Vec x
Vec b
KSP ksp
```

装配方式：

1. 为每个参与压力求解的流体粒子分配压力自由度。
2. 建立 `particle_index -> pressure_dof` 映射。
3. 预估每行非零元个数：
   - 普通流体行约为 `fluid_neighbor_count + 1`。
   - 自由面 Dirichlet 行为 1。
4. 创建 PETSc AIJ 矩阵并预分配。
5. 对每个流体粒子装配一行 PPE。
6. 调用 `MatAssemblyBegin/MatAssemblyEnd`。
7. 创建 RHS `Vec b` 和解向量 `Vec x`。
8. 使用 `KSP` 求解。

PETSc `MatSetValues` 使用 0-based 行列索引，符合当前 C++ 索引习惯。所有 `MatSetValues` 完成后必须执行矩阵 assembly。

## 5. PPE 行装配原则

### 5.1 自由面粒子

自由面压力 Dirichlet 条件直接替换方程行：

```text
A_ii = 1
b_i = P_fs
```

不使用 LSMPS pressure row。

### 5.2 非自由面流体粒子

PPE 原始形式：

```text
(1 / rho) laplacian(P) = (1 / dt) div(u*)
```

等价装配：

```text
laplacian(P) = rho / dt * div(u*)
```

对非自由面流体粒子，压力拉普拉斯使用 LSMPS pressure Neumann 逆矩阵：

```text
M = lsmps_matrices.particles[i].pressure_neumann.inverse_moment
```

流体邻居贡献：

```text
L_ij = (1 / re^2) * w_ij * (M_row_3 + M_row_4 + M_row_5) dot p_ij
```

矩阵项：

```text
A_i_i += -L_ij
A_i_j +=  L_ij
```

壁面 Neumann 源项：

```text
B_i = (1 / re) * sum_wall w_ij * (M_row_3 + M_row_4 + M_row_5) dot q_ij * C_j
```

其中：

```text
C_j = rho * g dot n_j
```

右端项：

```text
b_i = rho / dt * div(u*) - B_i
```

### 5.3 临时速度散度

PPE 右端项需要：

```text
div(u*) = dux/dx + duy/dy + duz/dz
```

第一版使用 LSMPS regular 逆矩阵计算：

```text
M = lsmps_matrices.particles[i].regular.inverse_moment
```

壁面粒子临时速度来自 Provisional 模块输出，不使用真实速度直接代替。

## 6. 模块边界

`pressure_poisson` 模块负责：

- 建立压力自由度映射。
- 使用 PETSc AIJ 构造 PPE 矩阵。
- 使用 PETSc Vec 构造 RHS。
- 调用 PETSc KSP 求解。
- 将 PETSc 解向量写回粒子压力数组或返回压力结果。
- 输出迭代次数、残差范数、收敛原因。

`pressure_poisson` 模块不负责：

- 计算临时速度。
- 更新粒子位置。
- 修正速度。
- 构造 LSMPS 逆矩阵。
- 自由面识别。

## 7. 建议接口

```cpp
struct PressurePoissonSolveInfo {
    bool converged = false;
    int iterations = 0;
    double final_residual_norm = 0.0;
    int converged_reason = 0;
};

struct PressurePoissonResult {
    std::vector<double> pressure;
    PressurePoissonSolveInfo solve_info;
};

class PressurePoissonAssembler {
public:
    PressurePoissonResult solve(
        const ParticleSet& particles,
        const TypedNeighborList& neighbors,
        const LsmpsMatrixSet& matrices,
        const ProvisionalVelocityResult& provisional,
        const SimulationConfig& config) const;
};
```

第一版可以只支持单进程 PETSc。接口和矩阵构造应避免阻碍后续 MPI 并行：

- 粒子自由度映射独立封装。
- 行装配逻辑独立封装。
- 矩阵预分配逻辑独立封装。

## 8. CMake 接入

建议通过 `pkg-config` 接入 PETSc：

```cmake
find_package(PkgConfig REQUIRED)
find_package(MPI REQUIRED)
pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc)

target_link_libraries(lsmps3d_core PUBLIC PkgConfig::PETSC MPI::MPI_CXX)
```

这样可以复用 Homebrew PETSc 提供的 include/lib 参数，同时补齐 MPI 头文件和链接参数。

本机 smoke 验证：

```text
mpicxx -std=c++17 petsc_smoke.cpp $(pkg-config --cflags --libs PETSc)
```

可以完成编译和运行。当前受限运行环境下 OpenMPI 可能打印 TCP bind 警告，但 PETSc 初始化和基本对象创建可用。

## 9. 第一版测试计划

### 9.1 PETSc smoke test

构造一个小型非对称稀疏线性系统，验证 PETSc 初始化、AIJ 矩阵装配、GMRES 求解和残差读取。

### 9.2 解析 Poisson 测试

构造规则粒子分布，给定解析压力场 `P_exact`，根据离散算子构造 RHS，再求解 PPE。验证：

- 压力误差。
- 残差。
- 迭代次数。
- Dirichlet 行是否正确。

### 9.3 静水压力测试

复用静水压力场：

```text
P = rho g (H - y)
```

检查 PPE 装配后的压力结果和理论静水压力是否一致。

### 9.4 VTK 输出

输出：

- `pressure`
- `pressure_error`
- `ppe_rhs`
- `ppe_residual`
- `ppe_converged_reason`
- `ppe_iterations`
- `lsmps_pressure_condition_number`

## 10. 后续优化

后续性能优化路径：

- PETSc MPI 多进程运行。
- 更准确的每行非零元预分配。
- 尝试 `fgmres`、`bcgs`、`lgmres`。
- 构建带 HYPRE 的 PETSc 并使用 BoomerAMG。
- 对自由面 Dirichlet 粒子采用已知值移项策略，减少耦合自由度。
- 分析矩阵条件数和粒子分布质量对收敛的影响。
