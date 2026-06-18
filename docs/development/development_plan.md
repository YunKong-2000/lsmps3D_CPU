# LSMPS 3D CPU 开发计划与进度记录

本文档用于记录三维自由面流动 LSMPS 程序的开发计划、当前完成情况、关键技术决策和后续待办事项。随着代码实现推进，应持续更新任务状态，避免设计文档和实际实现脱节。

## 1. 当前目标

当前阶段的目标是开发一个可在本地 CPU 上运行的三维 LSMPS 自由面流动求解器。第一版程序优先完成清晰、可验证的最小功能闭环：

1. 建立参数配置模块。
2. 生成三维粒子分布。
3. 输出 VTK 结果用于 ParaView 可视化。
4. 构建邻域搜索。
5. 识别自由面和近自由面粒子。
6. 构造 LSMPS type-A 基函数、核函数、常规逆矩阵和压力 Neumann 逆矩阵。
7. 完成逆矩阵、多项式重构和压力 Neumann 约束验证。
8. 计算临时速度 $\mathbf{u}^{*}$。
9. 组装并求解 PPE。
10. 使用压力梯度修正速度并更新粒子位置。

初始版本暂不追求复杂边界、并行性能和高级稳定化模型，重点是建立正确的程序结构和验证路径。

## 2. 已有理论与设计基础

当前已经具备以下文档基础：

- `docs/theory/lsmps3D.md`：程序总体框架、数据结构、模块划分、自由面和边界处理原则。
- `docs/theory/lsmps_scheme.md`：LSMPS type-A/type-B 数学格式、边界条件处理、PPE 离散与装配原则。

当前已明确的重要实现原则：

- PPE 构造不使用 LSMPS type-B 格式。
- PPE 压力拉普拉斯算子统一采用 type-A 格式。
- 自由面压力第一类边界条件不通过修改 LSMPS 算子实现。
- 自由面压力条件在压力线性系统中通过 Dirichlet 方程行直接施加。
- 壁面压力第二类边界条件可通过额外约束进入局部 LSMPS 矩矩阵。
- 自由面粒子属于流体粒子，不是与壁面粒子并列的顶层粒子类型。
- PPE 构造前必须先由 Provisional 模块计算临时速度 $\mathbf{u}^{*}$。
- Provisional 模块负责黏性项、重力项和外力项，不负责压力求解和最终位置更新。
- PPE 模块只负责根据 $\mathbf{u}^{*}$、自由面 Dirichlet 条件和壁面 Neumann 条件装配并求解压力。
- Correction 模块负责压力梯度、速度修正和下一时间步粒子位置更新。

## 3. 计划项目结构

初始代码结构建议如下：

```text
src/
    main.cpp
    core/
        vector3.*
        particle_set.*
        simulation_config.*
        simulation_state.*
    config/
        config_reader.*
        case_config.*
    geometry/
        domain_builder.*
        wall_builder.*
    neighbor/
        cell_linked_list.*
        neighbor_search.*
    free_surface/
        free_surface_detector.*
    lsmps/
        basis.*
        weight_function.*
        matrix_builder.*
        differential_operator.*
    provisional/
        provisional.*
    pressure_poisson/
        pressure_poisson.*
    correction/
        correction.*
    linear_solver/
        sparse_matrix.*
        cg_solver.*
    io/
        vtk_writer.*
        logger.*
tests/
    test_neighbor_search.*
    test_lsmps_operator.*
    test_provisional.*
    test_pressure_poisson.*
    test_correction.*
cases/
    hydrostatic_box/
    dam_break_3d/
```

实际开发中可以根据构建系统和依赖库选择调整文件名，但模块边界应尽量保持稳定。

当前项目骨架已按上述边界建立：

- `src/provisional/provisional.*`
- `src/pressure_poisson/pressure_poisson.*`
- `src/correction/correction.*`

## 4. 开发阶段

### 阶段 0：工程骨架

目标：建立可编译、可运行、可扩展的 C++ 工程基础。

状态：已完成

任务：

- [o] 确定构建系统，使用 CMake。
- [o] 建立 `src/`、`tests/`、`cases/`、`output/` 目录。
- [o] 建立基础 `main.cpp`。
- [o] 建立基础日志输出。
- [o] 确定 C++ 标准版本，当前使用 C++17。
- [o] 确定是否第一阶段引入 Eigen，当前阶段暂不引入。

验收标准：

- 可以完成一次空程序编译。
- 可以运行可执行文件并输出基础程序信息。
- 项目目录结构稳定。

### 阶段 1：核心数据结构

目标：实现粒子数据和仿真状态的基础表示。

状态：已完成

任务：

- [o] 实现三维向量或选择现有向量类型。
- [o] 实现 `ParticleSet`。
- [o] 实现顶层粒子类型枚举：`Fluid`、`Wall`。
- [o] 实现流体粒子状态枚举：`Internal`、`FreeSurface`、`NearFreeSurface`、`Splash`。
- [o] 实现 `SimulationState`。
- [o] 实现基础粒子初始化接口。

验收标准：

- 能够创建一组流体粒子和壁面粒子。
- 能够访问并修改位置、速度、压力、类型和流体状态。
- 数据结构能够支持后续邻域搜索和 VTK 输出。

### 阶段 2：参数配置模块

目标：优先建立统一的参数入口，避免后续模块依赖硬编码参数。

状态：已完成

任务：

- [o] 实现 `SimulationConfig`。
- [o] 确定第一阶段配置文件格式。
- [o] 支持读取时间步长、终止时间和输出间隔。
- [o] 支持读取粒子间距、支持域半径和初始几何参数。
- [o] 支持读取密度、黏度和重力。
- [o] 支持读取自由面识别相关阈值。
- [o] 支持读取线性求解器参数。
- [o] 提供默认参数和参数合法性检查。

当前实现说明：

- 第一阶段配置文件使用 INI 分块格式，支持 `#` 和 `;` 注释。
- 三维向量使用空格分隔，例如在 `[physical]` 中写 `gravity = 0.0 -9.81 0.0`。
- 程序启动支持 `lsmps3d [config_file]`；不提供配置文件时使用默认参数。
- 示例配置文件位于 `cases/hydrostatic_box/config.ini`。

验收标准：

- 程序启动时能够从配置文件或默认配置生成 `SimulationConfig`。
- 关键数值参数不再散落在代码中硬编码。
- 后续 VTK、邻域搜索、自由面识别和 LSMPS 算子都通过统一配置获取参数。

### 阶段 3：VTK 输出

目标：在进入邻域搜索前实现最小可用的 ParaView 可视化输出，尽早建立粒子数据检查手段。

状态：已完成

任务：

- [o] 输出粒子坐标。
- [o] 输出速度向量。
- [o] 输出压力。
- [o] 输出顶层粒子类型。
- [o] 输出流体粒子状态。
- [o] 输出邻居数量占位字段。
- [o] 输出算子或自由面识别诊断量的预留接口。

当前实现说明：

- VTK 输出采用 legacy ASCII `POLYDATA` 点云格式。
- 默认输出字段包括 `velocity`、`pressure`、`particle_type`、`fluid_state` 和 `neighbor_count`。
- `VtkWriter::writeParticles` 支持附加 scalar/vector 字段，供后续自由面识别和 LSMPS 诊断量复用。
- 程序启动后会写出一个最小检查文件 `output/initial_particles.vtk`。

验收标准：

- ParaView 能够打开输出文件。
- 可以按粒子类型、压力、速度和流体状态着色。
- 输出接口可被后续邻域搜索、自由面识别和算例复用，可以向vtk文件中灵活添加vector,scalar。

### 阶段 4：邻域搜索

目标：实现三维均匀网格邻域搜索。

状态：已完成

任务：

- [o] 实现 cell linked-list 数据结构。
- [o] 根据支持域半径 $r_e$ 构造背景网格。
- [o] 实现粒子到 cell 的映射。
- [o] 实现邻居查询。
- [o] 输出每个粒子的邻居数量。
- [o] 编写规则粒子分布下的邻居搜索测试。

当前实现说明：

- 邻域搜索采用哈希表形式的三维 cell linked-list，只存储非空 cell。
- cell 尺寸取支持域半径 `support_radius`，每个目标粒子检查自身 cell 及周围 26 个 cell。
- 搜索结果使用真实距离平方筛选，并包含距离正好等于支持半径的粒子。
- 流体粒子和壁面粒子统一参与邻域搜索。
- 已分别保存总邻居索引、流体邻居索引和壁面邻居索引。
- 已支持总邻居数、流体邻居数和壁面邻居数分类统计，并写入 VTK 默认字段。
- 已增加 50x50x50 规则流体粒子加单层壁面壳的大规模测试，支持域半径为 `3.1 * particle_spacing`，输出 `output/neighbor_search_50x50x50.vtk` 用于可视化验收。
- 实现细节见 `docs/development/modules/neighbor_search.md`。

验收标准：

- 对规则粒子分布，邻居数量和距离筛选结果正确。
- 必须支持流体粒子和壁面粒子共同参与邻域搜索。
- 邻域搜索模块不依赖具体物理方程。

### 阶段 5：自由面识别

目标：在 LSMPS 离散模块前完成自由面状态识别，为后续算子构造、PPE 边界条件和可视化提供粒子状态信息。

状态：进行中

算法设计：

- [o] 完成壁面感知的自由面识别算法设计，见 `docs/development/modules/free_surface_detection_design.md`。
- [o] 参考 `docs/development/modules/ALGORITHM_GPU_CN.md` 实现 `R_open + R_cone` 综合判定；当前版本不实现最大连通域判定。

任务：

- [o] 实现基于 cubed sphere 开口面积 `R_open` 和锥形开口 `R_cone` 的自由面识别。
- [o] 标记 `Internal`、`FreeSurface`、`NearFreeSurface`。
- [o] 设计并实现 `Splash` 状态的初始保守判据。
- [o] 将自由面识别阈值接入参数配置模块。
- [o] 将自由面状态和诊断量写入 VTK 输出用于测试验收。
- [o] 使用简单静态粒子分布检查自由面识别结果。

验收标准：

- 自由面粒子识别结果空间上连续且可视化合理。
- `FreeSurface`、`NearFreeSurface` 状态能够被后续 PPE 和诊断输出直接使用。
- 自由面识别模块只依赖粒子数据、邻居列表和配置参数，不依赖 LSMPS 离散算子。

### 阶段 6：LSMPS type-A 逆矩阵模块

目标：实现基础 LSMPS type-A 基函数、核函数、常规逆矩阵、压力 Neumann 逆矩阵和矩阵诊断接口，为后续物理离散模块提供公共能力。

状态：已完成

实现路径：

- [o] 建立 LSMPS type-A 逆矩阵模块实现路径文档，见 `docs/development/modules/lsmps_type_a_implementation.md`。

任务：

- [o] 实现二阶多项式基向量 $\mathbf{p}_{ij}$。
- [o] 实现壁面 Neumann 约束基向量 $\mathbf{q}_{ij}$。
- [o] 实现权函数 $w_{ij}$。
- [o] 引入线性代数库用于固定尺寸局部矩阵求逆和诊断。
- [o] 构造并缓存常规 `inverse_moment`，壁面邻居按 $\mathbf{p}_{ij}\mathbf{p}_{ij}^T$ 贡献。
- [o] 构造并缓存压力 Neumann `inverse_moment`，壁面邻居按 $\mathbf{q}_{ij}\mathbf{q}_{ij}^T$ 贡献。
- [o] 输出两类逆矩阵的秩、条件数、特征值和构造状态。
- [o] 对解析多项式验证常规逆矩阵重构能力。
- [o] 对解析压力 Neumann 约束验证压力逆矩阵构造能力。

当前实现说明：

- LSMPS 基函数、权函数和逆矩阵构造位于 `src/lsmps/`。
- 当前使用 Eigen 固定尺寸 9x9 矩阵完成特征值诊断、秩检查和求逆。
- 原始累加矩阵只在构造过程中临时存在，不长期缓存。
- 对每个流体粒子缓存 `regular.inverse_moment` 和 `pressure_neumann.inverse_moment`。
- `regular` 矩阵中壁面邻居按 $\mathbf{p}_{ij}\mathbf{p}_{ij}^T$ 贡献。
- `pressure_neumann` 矩阵中壁面邻居按 $\mathbf{q}_{ij}\mathbf{q}_{ij}^T$ 贡献。
- 测试 `lsmps_matrix_test` 会输出 `build/output/lsmps_matrix_diagnostics.vtk` 用于 ParaView 检查。
- 测试 `lsmps_operator_validation_test` 使用 20x20x20 流体粒子加单层壁面粒子的静水压力和管道解析速度算例，输出 `output/lsmps_hydrostatic_20x20x20.vtk` 与 `output/lsmps_pipe_20x20x20.vtk`。

验收标准：

- 对简单解析多项式，常规逆矩阵重构结果与解析多项式系数一致。
- 近壁流体粒子的压力 Neumann 逆矩阵能够稳定构造，并且不同于常规逆矩阵。
- 后续模块只需要读取缓存的 `inverse_moment`，不依赖长期缓存的原始矩阵。
- 壁面邻居默认参与矩阵构造，近壁粒子可输出壁面邻居数和两类矩阵质量诊断。
- 逆矩阵实现与 `docs/theory/lsmps_scheme.md` 中 type-A 及壁面第二类边界格式一致。

### 阶段 7：临时速度 Provisional 模块

目标：在 PPE 构造前，根据当前速度、黏性项、重力和外力计算临时速度 $\mathbf{u}^{*}$。

状态：未开始

任务：

- [ ] 设计 `src/provisional/provisional.*` 模块接口。
- [ ] 使用当前速度场计算黏性项。
- [ ] 加入重力项和后续可扩展的外力项。
- [ ] 计算流体粒子的临时速度 $\mathbf{u}^{*}$。
- [ ] 计算壁面粒子的临时速度，满足后续 PPE 右端项中壁面速度参与散度离散的需求。
- [ ] 输出临时速度、黏性加速度和外力加速度诊断字段。
- [ ] 使用解析速度场或简单黏性扩散算例验证临时速度更新。

验收标准：

- PPE 模块可以直接读取临时速度场计算 $\nabla \cdot \mathbf{u}^{*}$。
- 重力项和黏性项贡献可单独输出检查。
- 模块不负责压力求解，也不更新最终位置。

### 阶段 8：压力泊松方程

目标：实现 PPE 的矩阵组装和求解。

状态：未开始

任务：

- [ ] 设计稀疏矩阵数据结构或接入现有库。
- [ ] 实现 PPE 右端项 $\frac{\rho}{\Delta t}\nabla \cdot \mathbf{u}^{*}$。
- [ ] 使用 type-A 格式组装压力拉普拉斯算子。
- [ ] 实现自由面压力 Dirichlet 方程行。
- [ ] 实现壁面压力第二类边界贡献。
- [ ] 实现 CG 或其他线性求解器。
- [ ] 输出压力残差和迭代次数。
- [ ] 通过解析 Poisson 问题进行验证。

验收标准：

- PPE 不使用 type-B 格式。
- 自由面压力第一类边界通过线性系统行替换施加。
- 简单 Poisson 验证算例误差可控。

### 阶段 9：压力修正 Correction 模块

目标：在 PPE 得到压力后，计算压力梯度，对临时速度进行修正，并更新下一时间步粒子速度和位置。

状态：未开始

任务：

- [ ] 设计 `src/correction/correction.*` 模块接口。
- [ ] 使用 LSMPS 常规或压力相关逆矩阵计算压力梯度。
- [ ] 根据压力梯度修正临时速度，得到下一时间步速度。
- [ ] 根据修正后速度更新粒子位置。
- [ ] 处理自由面、近自由面、Splash 和壁面粒子在速度修正中的参与规则。
- [ ] 输出压力梯度、速度修正量、位移量和最大速度诊断字段。
- [ ] 记录 CFL、最大速度和压力残差。

验收标准：

- 临时速度、压力场和压力梯度能够闭合生成下一时间步速度。
- 粒子位置更新只发生在 correction 阶段。
- 模块不负责 PPE 线性系统装配和求解。

### 阶段 10：不可压缩流动时间推进闭环

目标：串联 Provisional、PPE 和 Correction，完成预测、压力求解、速度修正和位置更新的基本闭环。

状态：未开始

任务：

- [ ] 串联临时速度计算、PPE 装配求解和压力修正模块。
- [ ] 每个时间步更新邻域搜索、自由面识别和 LSMPS 逆矩阵。
- [ ] 管理时间步、输出间隔和诊断日志。
- [ ] 输出 CFL、最大速度、压力残差和自由面粒子数量。

验收标准：

- 程序能够完成多个时间步而不崩溃。
- 静水箱算例中压力分布接近静水压力。
- 输出文件可连续可视化。

### 阶段 11：基础自由面算例

目标：基于已有自由面识别、PPE 和时间推进模块完成第一个三维自由面流动算例。

状态：未开始

任务：

- [ ] 在 PPE 中对 `FreeSurface` 粒子施加压力 Dirichlet 条件。
- [ ] 完成三维水箱或三维溃坝初始算例。
- [ ] 输出自由面状态用于 ParaView 检查。

验收标准：

- 自由面粒子识别结果空间上连续且可视化合理。
- 自由面粒子的压力边界条件正确进入 PPE。
- 基础自由面算例能够稳定运行短时间。

## 5. 当前完成情况

| 日期 | 内容 | 状态 |
|---|---|---|
| 2026-06-17 | 建立程序总体框架文档 `docs/theory/lsmps3D.md` | 已完成 |
| 2026-06-17 | 整理 LSMPS 离散与 PPE 装配文档 `docs/theory/lsmps_scheme.md` | 已完成 |
| 2026-06-17 | 建立开发计划与进度记录文档 `docs/development/development_plan.md` | 已完成 |
| 2026-06-17 | 建立 CMake/C++17 工程骨架、基础日志模块、主程序和 smoke test | 已完成 |
| 2026-06-17 | 实现核心粒子数据结构、粒子类型枚举、流体状态枚举和基础测试 | 已完成 |
| 2026-06-18 | 实现参数配置模块、INI 配置读取、默认参数、合法性检查和配置测试 | 已完成 |
| 2026-06-18 | 实现 VTK legacy ASCII 粒子输出模块、默认字段和附加诊断字段接口 | 已完成 |
| 2026-06-18 | 建立 LSMPS type-A 逆矩阵模块实现路径文档 | 已完成 |
| 2026-06-18 | 实现 LSMPS type-A 逆矩阵模块、Eigen 接入、矩阵诊断和 VTK 测试输出 | 已完成 |
| 2026-06-18 | 增加 20x20x20 静水压力和管道解析速度算子验证，输出 VTK 可视化验收文件 | 已完成 |

## 6. 近期优先级

近期建议按以下顺序推进：

1. 建立 CMake 工程骨架。
2. 实现 `ParticleSet`、粒子类型枚举和基础状态数据。
3. 优先实现参数配置模块。
4. 实现 VTK 输出，尽早建立可视化反馈。
5. 实现邻域搜索。
6. 紧接邻域搜索实现自由面识别。
7. 开发 Provisional 模块，计算 PPE 所需临时速度。
8. 基于临时速度和 LSMPS 逆矩阵模块进入 PPE 组装和压力求解。
9. 开发 Correction 模块，使用压力梯度修正速度并更新粒子位置。
10. 串联 Provisional、PPE 和 Correction，形成不可压缩时间推进闭环。

这样可以尽早发现数据结构、参数传递、输出格式、邻域搜索和自由面判据中的问题，避免直接进入复杂流动求解后难以定位错误。

## 7. 待决策事项

- C++ 标准版本：建议 `C++17` 或更高。
- 线性代数库：工程骨架阶段暂不引入 Eigen；从 LSMPS 逆矩阵模块开始建议引入 Eigen，用于固定尺寸局部矩阵求逆、特征值诊断和后续 PPE 稀疏线性代数。
- 稀疏矩阵格式：初始阶段可使用 Eigen sparse 或自定义 CSR。
- VTK 格式：建议先支持 VTK legacy ASCII，后续再扩展二进制或 VTU。
- 配置文件格式：第一阶段采用 INI 分块格式，支持注释和三维向量；后续如算例复杂度提升，可再评估 JSON 或 YAML。
- 单元测试框架：可选择 Catch2、GoogleTest，或先用简单测试可执行文件。

## 8. 设计保留意见

### 8.1 粒子集合管理方式

当前实现中，流体粒子和壁面粒子统一存放在同一个 `ParticleSet` 中，并通过 `ParticleType` 区分 `Fluid` 和 `Wall`。这种方式便于早期实现邻域搜索、VTK 输出和统一索引管理，但仍存在需要后续复审的问题：

- 壁面粒子所需属性少于流体粒子，统一存储可能造成部分字段语义不清。
- 后续如果壁面粒子需要法向量、边界类型等专有属性，`ParticleSet` 可能变得臃肿。
- 如果流体粒子和壁面粒子在 PPE、LSMPS 矩阵构造或输出中长期采用不同访问模式，统一集合可能不利于模块边界清晰。

该问题暂不立即重构。建议在完成以下模块后重新评估：

1. VTK 输出模块。
2. 邻域搜索模块。
3. 自由面识别模块。
4. LSMPS type-A 算子初版。

复审时重点比较两种方案：

- 继续使用统一 `ParticleSet`，并通过更清晰的访问接口约束不同粒子类型可用字段。
- 拆分为流体粒子集合和壁面粒子集合，同时增加统一的空间点视图或邻域搜索索引映射。

## 9. 更新规则

每完成一个开发任务，应更新本文档：

- 将对应任务从 `[ ]` 改为 `[x]`。
- 在“当前完成情况”表中新增一条记录。
- 如果实现方案和文档不同，应同步更新 `docs/theory/lsmps3D.md` 或 `docs/theory/lsmps_scheme.md`。
- 如果发现新的数值稳定性问题，应在对应阶段下追加待办项。

本文档应作为开发过程中的主索引，记录“下一步做什么”和“已经完成了什么”。
