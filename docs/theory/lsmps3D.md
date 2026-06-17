# LSMPS 方法三维自由面流动模拟程序

## 1. 程序定位

本程序计划开发为一个基于 LSMPS（Least Squares Moving Particle Semi-implicit）方法的三维 CFD 求解器，主要面向不可压缩自由面流动问题。程序初期以 CPU 版本为主，优先保证算法结构清晰、模块边界明确和结果可验证，后续再逐步扩展并行计算、高阶模型和复杂边界能力。

目标应用包括：

- 三维水箱晃荡、溃坝、波浪传播等典型自由面问题。
- 含复杂固壁边界的不可压缩流动。
- 需要较强自由面变形能力的粒子法流动模拟。

初始版本建议以研究型求解器为定位，不急于追求完整工程软件能力。代码结构应便于检查每个数值环节的正确性。

## 2. 控制方程与基本假设

初始程序以三维不可压缩牛顿流体为主要对象，采用拉格朗日粒子描述。

连续方程：

$$
\nabla \cdot \mathbf{u} = 0
$$

动量方程：

$$
\frac{D\mathbf{u}}{Dt}
= -\frac{1}{\rho}\nabla p
+ \nu \nabla^2 \mathbf{u}
+ \mathbf{g}
+ \mathbf{f}_{\mathrm{ext}}
$$

其中：

- $\mathbf{u}$ 为速度矢量。
- $p$ 为压力。
- $\rho$ 为密度。
- $\nu$ 为运动黏度。
- $\mathbf{g}$ 为重力加速度。
- $\mathbf{f}_{\mathrm{ext}}$ 为可选外力项。

初始版本可采用如下假设：

- 单相不可压缩流体。
- 密度和黏度为常数。
- 自由面处压力取大气压，可设为零表压。
- 固壁采用单层壁面粒子表征，壁面粒子携带壁面法向向量，法向向量指向流体内部。
- 表面张力、湍流模型、可压缩效应暂不作为第一阶段目标。

## 3. 数值方法框架

### 3.1 LSMPS 空间离散

LSMPS 方法通过邻域粒子上的加权最小二乘重构计算空间导数。对目标粒子 $i$，在其支持域内选择邻居粒子 $j$，构造局部坐标差：

$$
\Delta x_{ij} = x_j - x_i,\qquad
\Delta y_{ij} = y_j - y_i,\qquad
\Delta z_{ij} = z_j - z_i
$$

对标量场 $\phi$，可在粒子 $i$ 附近进行局部泰勒展开：

$$
\phi_j - \phi_i
\approx
\frac{\partial \phi_i}{\partial x}\Delta x_{ij}
+ \frac{\partial \phi_i}{\partial y}\Delta y_{ij}
+ \frac{\partial \phi_i}{\partial z}\Delta z_{ij}
+ \frac{1}{2}\frac{\partial^2 \phi_i}{\partial x^2}\Delta x_{ij}^2
+ \frac{1}{2}\frac{\partial^2 \phi_i}{\partial y^2}\Delta y_{ij}^2
+ \frac{1}{2}\frac{\partial^2 \phi_i}{\partial z^2}\Delta z_{ij}^2
+ \cdots
$$

通过加权最小二乘问题求解局部导数：

$$
\min_{\mathbf{a}_i}
\sum_{j \in \mathcal{N}_i}
w_{ij}
\left(
\mathbf{p}_{ij}^{T}\mathbf{a}_i
- (\phi_j - \phi_i)
\right)^2
$$

其中，$\mathcal{N}_i$ 为粒子 $i$ 的邻居集合，$w_{ij}$ 为权函数，$\mathbf{p}_{ij}$ 为局部多项式基，$\mathbf{a}_i$ 为待求导数系数。

对标量场或矢量场进行局部多项式拟合，从而得到：

- 一阶导数：$\partial f / \partial x$，$\partial f / \partial y$，$\partial f / \partial z$
- 二阶导数：$\partial^2 f / \partial x^2$，$\partial^2 f / \partial y^2$，$\partial^2 f / \partial z^2$
- 混合导数：视精度需求决定是否保留
- 梯度、散度、拉普拉斯算子

第一阶段建议实现二阶精度的 LSMPS 算子，不需要封装成固定接口：

$$
\nabla \phi,\qquad
\nabla \cdot \mathbf{u},\qquad
\nabla^2 \phi
$$

为便于调试，应保留每个粒子的矩阵条件数、邻居数量和算子构造状态，用于检查粒子分布质量。

### 3.2 邻域搜索

三维粒子计算中邻域搜索是核心性能瓶颈。初始版本建议使用均匀网格或 cell linked-list 方法：

- 根据支持域半径 $r_e$ 划分背景网格。
- 每个时间步更新粒子所在 cell。
- 搜索当前 cell 及周围相邻 cell 内的候选粒子。
- 根据距离筛选有效邻居。

邻域搜索模块应独立于具体物理模型，供 LSMPS 算子、自由面识别、边界处理共同使用。

### 3.3 时间推进

初始版本可采用半隐式投影法结构。一个典型时间步如下：

1. 更新邻域搜索结构。
2. 识别自由面粒子和近自由面粒子。
3. 构造或更新 LSMPS 微分使用到的矩阵。
4. 计算黏性项、重力项和外力项。
5. 预测中间速度 $\mathbf{u}^{*}$。
6. 根据 $\mathbf{u}^{*}$ 和不可压缩条件构造压力泊松方程。
7. 求解压力 $p$。
8. 使用压力梯度修正速度。
9. 更新粒子位置。
10. 处理自由面、固壁边界和粒子分布修正。
11. 输出结果和诊断信息。

可写成伪代码：

```text
initialize_case()

for step in time_steps:
    build_neighbor_list()
    build_lsmps_operators()

    apply_boundary_condition()
    detect_free_surface()

    compute_explicit_terms()
    predict_velocity()

    assemble_pressure_poisson()
    solve_pressure()

    correct_velocity()
    update_particle_position()

    stabilize_particle_distribution()
    write_output()
```

## 4. 自由面处理

自由面流动是本程序的主要目标之一。初始版本需要至少包含自由面粒子识别和自由面压力边界条件。需要注意，自由面粒子仍然是流体粒子，不应作为与壁面粒子并列的独立粒子类型；它只是流体粒子在当前时间步的一种状态。

### 4.1 自由面识别

可考虑以下判据：

- 邻居数量低于内部粒子阈值。
- 粒子数密度低于参考值。
- LSMPS 局部矩阵质量异常。
- 法向量或局部空缺方向判据。

第一阶段可采用简单稳定的邻居数量或粒子数密度判据，后续再引入更鲁棒的几何识别方法。识别结果建议写入流体粒子状态枚举，至少包含内部粒子、自由面粒子、近自由面粒子和飞溅粒子。

### 4.2 自由面压力条件

对被识别为自由面状态的流体粒子，可设置：

$$
p = 0
$$

或在压力泊松方程中作为 Dirichlet 条件处理。需要注意自由面识别的抖动可能导致压力场噪声，因此应保留平滑或滞回机制作为后续扩展点。

### 4.3 粒子分布修正

自由面和强变形区域容易出现粒子聚集或空洞。可在后续阶段加入：

- 粒子位移修正（particle shifting）。
- 局部重采样。
- 基于数密度的弱修正。
- 自由面附近的特殊限制，避免破坏自由面形状。

初始版本可先实现最小修正策略，并通过诊断量监控粒子质量。

## 5. 边界条件

程序需要支持以下边界类型：

- 固壁边界：无滑移。
- 自由面边界：压力 Dirichlet 条件。
- 开边界：作为后续扩展。
- 运动边界：作为后续扩展。

固壁处理可以从壁面粒子方法开始：

- 壁面粒子参与邻域搜索和 LSMPS 算子构造。

## 6. 压力泊松方程

不可压缩条件通过压力泊松方程实现。典型形式为：

$$
\nabla^2 p
=
\frac{\rho}{\Delta t}
\nabla \cdot \mathbf{u}^{*}
$$

压力方程的离散由 LSMPS 拉普拉斯算子给出。初始版本建议采用稀疏矩阵组装方式：

- 每个参与压力求解的流体粒子对应一个压力未知量。
- 内部粒子和近自由面粒子使用泊松方程。
- 自由面流体粒子施加压力 Dirichlet 条件。
- 固壁边界通过 Neumann 条件处理。

线性求解器可从共轭梯度法、BiCGSTAB 或 Eigen/自写稀疏求解器开始。后续如需提升性能，可接入更成熟的线性代数库。

## 7. 数据结构设计

建议初始阶段以结构清晰为优先，核心数据结构包括：

### 7.1 粒子数据

```text
ParticleSet
    position[N][3]
    velocity[N][3]
    pressure[N]
    density[N]
    type[N]
    fluid_state[N]
    mass[N]
    volume[N]
    neighbor_list[N]
```

顶层粒子类型用于区分物理角色，建议至少包含：

```text
Fluid
Wall
```

其中，`Fluid` 粒子需要进一步保存流体粒子状态枚举：

```text
Internal
FreeSurface
NearFreeSurface
Splash
```

含义如下：

- `Internal`：内部流体粒子，按完整不可压缩流体方程参与计算。
- `FreeSurface`：自由面流体粒子，通常施加自由面压力条件。
- `NearFreeSurface`：近自由面流体粒子，用于增强自由面附近判据和数值处理的连续性。
- `Splash`：飞溅流体粒子，可在压力求解、邻域构造或时间步限制中采用特殊处理。

如果使用 C++ 实现，建议优先采用结构化数组（SoA）布局，便于后续优化缓存访问和并行化。

### 7.2 计算配置

```text
SimulationConfig
    dt
    end_time
    output_interval
    particle_spacing
    support_radius
    density
    viscosity
    gravity
    lsmps_order
    case_name
```

### 7.3 求解状态

```text
SimulationState
    current_step
    current_time
    max_velocity
    cfl_number
    pressure_residual
    neighbor_statistics
```

这些状态量用于时间步控制、日志输出和结果诊断。

## 8. 模块划分

建议初始代码按以下模块组织：

```text
src/
    main.cpp
    core/
        particle_set.*
        simulation_config.*
        simulation_state.*
    geometry/
        domain_builder.*
        boundary_geometry.*
    neighbor/
        cell_linked_list.*
        neighbor_search.*
    lsmps/
        basis.*
        weight_function.*
        operator_builder.*
        differential_operator.*
    physics/
        incompressible_solver.*
        pressure_poisson.*
        free_surface.*
        boundary_condition.*
        particle_shifting.*
    linear_solver/
        sparse_matrix.*
        cg_solver.*
        bicgstab_solver.*
    io/
        case_reader.*
        vtk_writer.*
        logger.*
    tests/
        test_lsmps_operator.*
        test_neighbor_search.*
        test_pressure_solver.*
```

模块职责：

- `core`：基础数据结构，不包含具体数值算法。
- `geometry`：初始粒子生成、边界几何和算例定义。
- `neighbor`：邻居搜索和空间哈希。
- `lsmps`：LSMPS 矩阵构造、权函数和微分算子。
- `physics`：不可压缩流体时间推进、自由面和边界处理。
- `linear_solver`：压力方程相关线性求解器。
- `io`：输入参数、日志和 VTK/CSV 输出。
- `tests`：单元测试和基础验证算例。

## 9. 输入输出

### 9.1 输入

初始阶段可以使用简单文本、JSON 或 YAML 配置文件描述算例：

```text
case_name
domain_size
particle_spacing
time_step
end_time
fluid_properties
boundary_type
output_interval
```

几何生成可先从内置算例开始，例如：

- 静水箱。
- 三维溃坝。
- 水柱坍塌。
- 简单波浪传播。

### 9.2 输出

建议优先支持 VTK 格式，便于使用 ParaView 后处理。输出字段包括：

- 粒子位置。
- 速度。
- 压力。
- 粒子类型。
- 流体粒子状态。
- 邻居数量。
- 数值诊断量。

日志应记录：

- 当前步数和物理时间。
- 最大速度。
- CFL 数。
- 压力迭代步数和残差。
- 自由面粒子数量。
- 最小、最大、平均邻居数量。

## 10. 验证算例

为了逐步验证程序，建议按以下顺序推进：

1. 邻域搜索验证：规则粒子分布下邻居数量是否正确。
2. LSMPS 导数验证：对解析函数计算梯度和拉普拉斯误差。
3. 压力泊松验证：给定解析解的 Poisson 方程。
4. 静水压力验证：三维水箱内压力随深度线性变化。
5. 溃坝算例：比较自由面形态和前沿位置。
6. 水箱晃荡：验证自由面大变形和壁面压力。

每个验证算例都应保存输入参数、输出指标和误差评估脚本，避免后续修改破坏已有能力。

## 11. 开发阶段规划

### 阶段一：最小可运行求解器

- 建立项目结构和构建系统。
- 实现粒子数据结构。
- 实现三维均匀网格邻域搜索。
- 实现 LSMPS 梯度和拉普拉斯算子。
- 输出 VTK 文件。
- 完成解析函数导数验证。

### 阶段二：不可压缩流动核心

- 实现半隐式时间推进。
- 组装压力泊松方程。
- 实现基础线性求解器。
- 实现速度修正。
- 完成静水压力验证。

### 阶段三：自由面流动

- 实现自由面粒子识别。
- 实现自由面压力边界条件。
- 实现基础粒子分布修正。
- 完成三维溃坝和水柱坍塌算例。

### 阶段四：边界和稳定性增强

- 改进固壁边界处理。
- 加入更稳健的自由面识别。
- 增加自适应时间步。
- 增加诊断工具和错误检查。

### 阶段五：性能优化和扩展

- 优化数据布局和缓存访问。
- 加入 OpenMP 多线程。
- 优化稀疏矩阵组装和线性求解。
- 支持更复杂几何输入。

## 12. 初始实现原则

为了降低开发风险，初始版本建议遵循以下原则：

- 先保证三维框架正确，再优化性能。
- 每个数值算子都必须有独立测试。
- 自由面、边界、压力求解应尽量解耦。
- 不在早期引入过多经验修正项。
- 所有关键参数都应能从配置文件读取。
- 输出足够多的诊断量，便于定位数值不稳定来源。
- 代码结构应为后续并行化预留空间。

## 13. 待进一步明确的问题

后续需要在本文档基础上继续细化的问题包括：

- 程序主要使用 C++、Fortran 还是其他语言实现：程序使用c++实现在本地运行。
- 是否依赖 Eigen、PETSc、Trilinos 等外部线性代数库：尽可能使用成熟的现有数学计算库，避免重复造轮子。
- LSMPS 使用的多项式阶数、权函数形式和支持域半径：会提供额外的数学格式文档。
- 压力泊松方程的边界条件具体离散方式：使用单层壁面粒子表征壁面，同上，边界条件施加方法会在数学格式文档中给出。
- 自由面识别采用单一判据还是组合判据：尽可能简化和高效，减少可调参数。
- 固壁边界采用壁面粒子、虚粒子还是镜像粒子：单层壁面粒子
- 是否需要从第一版开始支持 OpenMP：暂时不需要，视情况而定。
- 结果输出使用 VTK legacy、VTU 还是 HDF5/XDMF：一定要支持vtk格式便于在paraview中可视化计算结果。

本文档作为初始框架，后续应随着算法选择和代码实现逐步更新为更具体的设计说明。
