# LSMPS Type-A 逆矩阵模块实现路径

本文档描述 LSMPS type-A 基础模块的第一版实现路径。该模块负责提供基函数、核函数、局部矩阵求逆结果和矩阵质量诊断，为后续速度散度、压力拉普拉斯、黏性项、压力梯度等模块复用。

当前状态：第一版已实现。

代码位置：

- `src/lsmps/lsmps_basis.*`
- `src/lsmps/weight_function.*`
- `src/lsmps/lsmps_matrices.*`
- `tests/lsmps_matrix_test.cpp`
- `tests/lsmps_operator_validation_test.cpp`

测试输出：

- `build/output/lsmps_matrix_diagnostics.vtk`
- `output/lsmps_hydrostatic_20x20x20.vtk`
- `output/lsmps_pipe_20x20x20.vtk`

为避免符号混乱，本文档采用以下实现命名：

- `raw_moment`：原始累加矩阵，例如 `sum(w p p^T)` 或 `sum(w p p^T) + sum(w q q^T)`。
- `inverse_moment`：`raw_moment^{-1}`。这对应理论文档中实际进入离散算子各行的矩阵。

后续离散算子使用的是 `inverse_moment` 的行向量，不需要长期访问 `raw_moment`。因此第一版实现只缓存 `inverse_moment` 和诊断量；`raw_moment` 只在构造、求逆和诊断过程中临时存在。

核心原则：

- 壁面粒子默认参与 LSMPS 相关矩阵构造。
- 普通标量/速度类算子使用常规 type-A 逆矩阵，壁面邻居按普通几何邻居贡献 `p p^T`。
- 压力相关算子需要代入壁面压力第二类边界条件，因此必须构造单独的压力 Neumann 逆矩阵；壁面邻居贡献 `q q^T`，不是 `p p^T`。
- 压力 Neumann 逆矩阵在 LSMPS 模块中统一计算和缓存，不在 PPE 模块中临时现算。
- LSMPS 模块不负责装配 PPE 稀疏矩阵，但应提供 PPE 装配所需的压力专用逆矩阵和壁面 `q` 基函数。
- 当前程序不采用 type-B 构造 PPE。

## 1. 模块目标

阶段 6 的目标是实现一个可复用、可诊断、GPU 可迁移的 LSMPS type-A 逆矩阵模块：

- 二阶 3D type-A 基向量 `p_ij`。
- 壁面 Neumann 约束基向量 `q_ij`。
- 紧支撑核函数/权函数 `w_ij`。
- 常规逆矩阵 `regular.inverse_moment`。
- 压力 Neumann 逆矩阵 `pressure_neumann.inverse_moment`。
- 每类逆矩阵对应的状态、邻居数、秩、特征值和条件数诊断。
- 解析多项式重构测试和近壁压力矩阵测试。

本模块第一版不直接对外提供固定的梯度、散度、拉普拉斯函数。具体算子的系数提取和物理右端项仍放在后续模块中实现，但这些模块必须复用本模块已经缓存好的逆矩阵。

## 2. 矩阵类型

### 2.1 常规逆矩阵

常规矩阵用于普通标量场、速度散度、黏性项等不需要压力壁面第二类边界约束的算子。

对目标流体粒子 `i`：

```text
raw_regular_i =
    sum_{j in Fluid, j != i} w_ij p_ij p_ij^T
  + sum_{j in Wall}          w_ij p_ij p_ij^T

regular_inverse_i = raw_regular_i^{-1}
```

这里壁面粒子按普通几何邻居参与。对于速度散度，壁面粒子的临时速度如何取值由后续时间推进或边界模块决定，LSMPS 模块只负责几何逆矩阵。

### 2.2 压力 Neumann 逆矩阵

压力相关算子在近壁区域需要代入壁面压力第二类边界条件：

```text
grad(P) dot n = C
```

因此其局部最小二乘矩阵与常规矩阵不同。对目标流体粒子 `i`：

```text
raw_pressure_i =
    sum_{j in Fluid, j != i} w_ij p_ij p_ij^T
  + sum_{j in Wall}          w_ij q_ij q_ij^T

pressure_neumann_inverse_i = raw_pressure_i^{-1}
```

其中 `q_ij` 来自壁面法向和 type-A 基函数梯度：

```text
q_ij = r_e grad(p_ij) dot n_j
```

压力 PPE 模块后续使用 `pressure_neumann_inverse_i` 装配压力拉普拉斯系数和壁面边界源项，不再重算该逆矩阵。

当粒子没有壁面邻居时，`raw_pressure_i` 退化为只包含流体邻居的压力矩阵。实现上仍可以统一生成 `pressure_neumann_inverse_i`，这样 PPE 模块只需读取同一缓存。

## 3. 数学格式

对目标粒子 `i` 和邻居粒子 `j`：

```text
dx = x_j - x_i
dy = y_j - y_i
dz = z_j - z_i
re = support_radius
```

二阶 3D type-A 基向量长度为 9：

```text
p = [
    dx / re,
    dy / re,
    dz / re,
    dx^2 / (2 re^2),
    dy^2 / (2 re^2),
    dz^2 / (2 re^2),
    dx dy / re^2,
    dx dz / re^2,
    dy dz / re^2
]
```

壁面 Neumann 约束基向量长度同样为 9：

```text
q = [
    n_x,
    n_y,
    n_z,
    n_x dx / re,
    n_y dy / re,
    n_z dz / re,
    (n_x dy + n_y dx) / re,
    (n_x dz + n_z dx) / re,
    (n_y dz + n_z dy) / re
]
```

这里 `n` 是壁面粒子法向，方向约定必须与理论文档和边界模块一致。

如果某个后续模块需要从导数特征向量 `D_i` 中读取导数，分量约定为：

```text
D = [
    re * dphi/dx,
    re * dphi/dy,
    re * dphi/dz,
    re^2 * d2phi/dx2,
    re^2 * d2phi/dy2,
    re^2 * d2phi/dz2,
    re^2 * d2phi/dxdy,
    re^2 * d2phi/dxdz,
    re^2 * d2phi/dydz
]
```

## 4. 权函数

第一版使用紧支撑、实现简单的线性权函数：

```text
q = r / re
w(q) = max(0, 1 - q)
```

其中：

```text
r = sqrt(dx^2 + dy^2 + dz^2)
```

若 `r <= eps` 或 `r > re`，跳过该邻居。

后续可增加可配置核函数，例如 `linear`、`quadratic`、`cubic spline` 或 `Wendland`。

## 5. 数据结构设计

第一版建议引入 Eigen 作为固定尺寸 9x9 局部矩阵的线性代数后端。

建议类型：

```cpp
constexpr int basis_size = 9;

using BasisVector = Eigen::Matrix<double, basis_size, 1>;
using MomentMatrix = Eigen::Matrix<double, basis_size, basis_size>;
```

每一类逆矩阵的缓存：

```cpp
enum class LsmpsMatrixStatus {
    Valid,
    NotEnoughNeighbors,
    RankDeficient,
    IllConditioned,
    InversionFailed,
};

struct LsmpsInverseMatrix {
    MomentMatrix inverse_moment;
    LsmpsMatrixStatus status;
    int rank;
    double condition_number;
    double min_eigenvalue;
    double max_eigenvalue;
    std::size_t total_neighbor_count;
    std::size_t fluid_neighbor_count;
    std::size_t wall_neighbor_count;
};

struct LsmpsParticleMatrices {
    LsmpsInverseMatrix regular;
    LsmpsInverseMatrix pressure_neumann;
};

struct LsmpsMatrixSet {
    std::vector<LsmpsParticleMatrices> particles;
};
```

说明：

- 不缓存 `raw_moment`。
- `inverse_moment` 是后续离散算子真正需要的矩阵。
- `regular` 和 `pressure_neumann` 必须分开存储，因为壁面邻居的矩阵贡献不同。
- 诊断量来自临时 `raw_moment` 的谱分析，计算完成后保留数值诊断即可。

## 6. 为什么使用线性代数库

LSMPS 模块仍然需要构造 9x9 临时矩阵并求逆。手写 Gauss-Jordan 的依赖少，但不适合作为长期实现：

- 数值稳定性、主元策略和异常处理都需要自行维护。
- 很难可靠地区分奇异、近奇异、秩亏和病态矩阵。
- 后续 PPE 阶段仍然需要稀疏线性代数能力。
- 手写求逆会让维护重点偏离 LSMPS 算法本身。

使用 Eigen 或类似线性代数库的收益：

- 固定尺寸 9x9 矩阵开销可控。
- 可用 `SelfAdjointEigenSolver` 获取特征值并计算条件数。
- 可用 `FullPivLU` 或 `CompleteOrthogonalDecomposition` 判断秩。
- 可通过分解求解单位矩阵得到 `inverse_moment`，避免自己维护求逆流程。

线性代数库不会自动判断“物理坏粒子”，但能提供可靠的数值诊断。推荐判据：

```text
neighbor_count < min_neighbors    -> NotEnoughNeighbors
rank < 9                         -> RankDeficient
min_eigenvalue <= eigen_tol       -> RankDeficient 或 IllConditioned
condition_number > fail_threshold -> InversionFailed 或 IllConditioned
condition_number > warn_threshold -> IllConditioned
```

第一版可以对 `IllConditioned` 仍缓存逆矩阵，但必须在 VTK 中输出状态，便于检查坏粒子区域。

## 7. 构造流程

对每个流体目标粒子 `i`，分别构造两套临时矩阵。

### 7.1 常规逆矩阵

```text
raw = zero 9x9

for j in neighbors.fluid[i]:
    p = evaluateTypeABasis(x_j - x_i, re)
    raw += w_ij p p^T

for j in neighbors.wall[i]:
    p = evaluateTypeABasis(x_j - x_i, re)
    raw += w_ij p p^T

regular.inverse_moment = inverse(raw)
regular.diagnostics = diagnose(raw)
discard raw
```

### 7.2 压力 Neumann 逆矩阵

```text
raw = zero 9x9

for j in neighbors.fluid[i]:
    p = evaluateTypeABasis(x_j - x_i, re)
    raw += w_ij p p^T

for j in neighbors.wall[i]:
    q = evaluateTypeANeumannBasis(x_j - x_i, wall_normal_j, re)
    raw += w_ij q q^T

pressure_neumann.inverse_moment = inverse(raw)
pressure_neumann.diagnostics = diagnose(raw)
discard raw
```

注意：

- 两类矩阵都默认考虑壁面邻居，但壁面贡献形式不同。
- 常规矩阵的壁面贡献是 `p p^T`。
- 压力 Neumann 矩阵的壁面贡献是 `q q^T`。
- 最小邻居数和秩诊断应分别针对两类矩阵独立计算。
- 自由面压力 Dirichlet 条件不进入这两类矩阵，仍在 PPE 系统装配时替换方程行。

## 8. 提供给后续模块的接口

LSMPS 模块应提供逆矩阵缓存和低层基函数接口：

```cpp
BasisVector evaluateTypeABasis(const Vector3& offset, double support_radius);

BasisVector evaluateTypeANeumannBasis(
    const Vector3& offset,
    const Vector3& wall_normal,
    double support_radius);

double evaluateWeight(double distance, double support_radius, KernelType kernel);

LsmpsMatrixSet buildLsmpsMatrices(
    const ParticleSet& particles,
    const TypedNeighborList& neighbors,
    const LsmpsConfig& config);
```

后续普通算子模块读取：

```text
matrices.particles[i].regular.inverse_moment
```

后续压力相关模块读取：

```text
matrices.particles[i].pressure_neumann.inverse_moment
```

PPE 模块仍然负责：

- 提取压力拉普拉斯行系数。
- 计算壁面 Neumann 源项中的 `C_j`。
- 装配稀疏线性系统。
- 对自由面粒子施加 Dirichlet 方程行。

PPE 模块不负责：

- 重建压力 Neumann `raw_moment`。
- 求逆压力 Neumann 矩阵。
- 重新诊断局部 LSMPS 矩阵质量。

## 9. GPU 可迁移性

当前 CPU 实现可以使用 Eigen，但算法组织应保持 GPU 可迁移：

### Pass 1：构造并求逆常规矩阵

```text
one thread / one small thread group per particle
read fluid neighbor range
read wall neighbor range
accumulate regular raw_moment
diagnose and invert
write regular inverse_moment and diagnostics
```

### Pass 2：构造并求逆压力 Neumann 矩阵

```text
one thread / one small thread group per particle
read fluid neighbor range
read wall neighbor range
accumulate pressure raw_moment with p for fluid and q for wall
diagnose and invert
write pressure_neumann inverse_moment and diagnostics
```

GPU 版本应使用分开的 CSR 邻居列表：

```text
fluid_neighbor_offsets[N + 1]
fluid_neighbor_indices[Mf]
wall_neighbor_offsets[N + 1]
wall_neighbor_indices[Mw]
```

这和当前 `TypedNeighborList` 中流体邻居、壁面邻居分开存储的方向一致。

## 10. VTK 诊断输出

LSMPS 模块应分别输出常规矩阵和压力 Neumann 矩阵的诊断字段：

```text
lsmps_regular_status
lsmps_regular_rank
lsmps_regular_condition_number
lsmps_regular_min_eigenvalue
lsmps_regular_max_eigenvalue

lsmps_pressure_status
lsmps_pressure_rank
lsmps_pressure_condition_number
lsmps_pressure_min_eigenvalue
lsmps_pressure_max_eigenvalue

lsmps_total_neighbor_count
lsmps_fluid_neighbor_count
lsmps_wall_neighbor_count
```

这些字段用于 ParaView 检查：

- 自由面附近矩阵是否退化。
- 近壁区域常规矩阵和压力矩阵的质量差异。
- 压力 Neumann 矩阵是否因为 `q q^T` 约束获得足够秩。
- 稀疏区域、壁面交界和角点处是否出现病态矩阵。

## 11. 测试计划

### 11.1 基础矩阵测试

构造规则 3D 粒子块，使用足够大的支持半径。

检查：

- 内部流体粒子的 `regular` 和 `pressure_neumann` 状态为 `Valid`。
- 无壁面邻居的内部粒子，两类逆矩阵应一致或近似一致。
- 无崩溃、无 NaN。

### 11.2 近壁压力矩阵测试

构造带单层壁面粒子的规则水箱。

检查近壁流体粒子：

- `wall_neighbor_count > 0`。
- `regular.inverse_moment` 与 `pressure_neumann.inverse_moment` 不应完全相同。
- `pressure_neumann` 的秩、条件数、最小特征值输出合理。
- VTK 中可以区分两类矩阵状态。

### 11.3 多项式重构测试

用常规逆矩阵验证普通 type-A 多项式重构：

```text
phi(x, y, z) = 2x - 3y + 4z + 5
phi(x, y, z) = x^2 + 2y^2 + 3z^2 + xy - 2xz + yz
```

测试流程：

```text
rhs = sum(w p (phi_j - phi_i))
D = regular.inverse_moment * rhs
compare D with analytical coefficients
```

### 11.4 压力 Neumann 约束测试

构造简单壁面法向和解析压力场，使 `grad(P) dot n = C` 可解析。

测试流程：

```text
rhs =
    sum_{fluid} w p (P_j - P_i)
  + sum_{wall}  w q re C_j

D = pressure_neumann.inverse_moment * rhs
compare normal derivative and selected derivatives with analytical values
```

### 11.5 坏粒子诊断测试

构造退化邻域，例如：

- 邻居全部近似共面。
- 邻居数量不足。
- 邻居分布极端偏向一侧。
- 壁面法向异常或角点区域。

检查 `NotEnoughNeighbors`、`RankDeficient`、`IllConditioned` 是否能按预期出现，并输出到 VTK 中可视化。

### 11.6 20x20x20 静水压力算子验收

构造 20x20x20 规则流体粒子块，并在外围添加单层壁面粒子。给定理论静水压力：

```text
P = rho g (H - y)
```

使用 `pressure_neumann.inverse_moment` 计算压力梯度和压力拉普拉斯。理论结果为：

```text
grad(P) = (0, -rho g, 0)
laplacian(P) = 0
```

输出文件：

```text
output/lsmps_hydrostatic_20x20x20.vtk
```

主要可视化字段：

- `pressure_gradient`
- `pressure_gradient_error`
- `pressure_laplacian`
- `pressure_laplacian_error`
- `lsmps_pressure_status`
- `lsmps_pressure_rank`
- `lsmps_pressure_condition_number`
- `lsmps_wall_neighbor_count`

### 11.7 20x20x20 管道解析速度算子验收

构造 20x20x20 规则流体粒子块，并在外围添加单层壁面粒子。给定无散解析速度场：

```text
u_x = y^2 - z^2
u_y = 2 z^2 - x^2
u_z = 3 x^2 - y^2
```

理论结果为：

```text
div(u) = 0
laplacian(u) = (0, 2, 4)
```

使用 `regular.inverse_moment` 计算速度梯度、速度散度和速度拉普拉斯。

输出文件：

```text
output/lsmps_pipe_20x20x20.vtk
```

主要可视化字段：

- `velocity_grad_xx` 到 `velocity_grad_zz`
- `velocity_divergence`
- `velocity_divergence_error`
- `velocity_laplacian`
- `velocity_laplacian_x`
- `velocity_laplacian_y`
- `velocity_laplacian_z`
- `velocity_laplacian_error`
- `lsmps_regular_status`
- `lsmps_regular_rank`
- `lsmps_regular_condition_number`
- `lsmps_wall_neighbor_count`

## 12. 第一版实现顺序

建议按以下顺序开发：

1. 引入 Eigen 依赖。
   - CMake 中使用 `find_package(Eigen3 REQUIRED)`。
   - `lsmps3d_core` 链接 `Eigen3::Eigen`。
2. `lsmps_basis.*`
   - `evaluateTypeABasis(offset, support_radius)`。
   - `evaluateTypeANeumannBasis(offset, wall_normal, support_radius)`。
3. `weight_function.*`
   - `evaluateWeight(distance, support_radius, kernel_type)`。
4. `moment_matrix.*`
   - 临时 9x9 `raw_moment` 构造。
   - 特征值、秩和条件数诊断。
   - 通过线性代数库求 `inverse_moment`。
   - 不缓存 `raw_moment`。
5. `lsmps_matrix_builder.*`
   - 构造 `regular.inverse_moment`。
   - 构造 `pressure_neumann.inverse_moment`。
   - 输出两类矩阵诊断。
6. `lsmps_matrix_test.cpp`
   - 常规矩阵构造。
   - 近壁压力 Neumann 矩阵构造。
   - 多项式重构。
   - 压力 Neumann 解析约束。
   - 坏粒子诊断。
   - VTK 诊断输出。

## 13. 当前不做的内容

第一版明确不做：

- 不实现 type-B。
- 不在 LSMPS 模块中装配 PPE 稀疏矩阵。
- 不在 LSMPS 模块中施加自由面压力 Dirichlet 方程行。
- 不把梯度、散度、拉普拉斯作为本模块的正式对外算子。
- 不长期缓存 `raw_moment`。
- 不手写 9x9 稠密矩阵求逆作为默认方案。
- 不实现 GPU 代码，只保持实现路径可迁移。

这些内容分别属于后续物理离散模块、PPE 装配阶段或 GPU 化阶段。
