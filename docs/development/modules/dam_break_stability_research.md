# 溃坝算例稳定性问题与改进方向研究

本文档记录 3D 溃坝算例在未使用 particle shifting 时出现的前沿粒子分布退化问题，并总结后续需要开发的稳定性改进方向。

## 1. 当前观测现象

算例输出目录：

- `build/output/dam_break_3d`

当前输出文件采用流体和壁面分离方式：

- `dam_break_3d_fluid_00000.vtk` 至 `dam_break_3d_fluid_00022.vtk`
- `dam_break_3d_wall_00000.vtk`

流体 VTK 文件当前包含字段：

- `POINTS`：粒子位置
- `velocity`
- `pressure`
- `fluid_state`
- `neighbor_count`
- `pressure_gradient`

阶段一诊断增强后，流体 VTK 还包含：

- `nearest_fluid_distance`
- `max_fluid_neighbor_distance`
- `mean_fluid_neighbor_distance`
- `number_density`
- `number_density_ratio`
- `x_positive_neighbor_count`
- `x_negative_neighbor_count`
- `y_positive_neighbor_count`
- `y_negative_neighbor_count`
- `z_positive_neighbor_count`
- `z_negative_neighbor_count`
- `directional_coverage_score`
- `geometry_min_eigenvalue`
- `geometry_max_eigenvalue`
- `geometry_condition_number`

从 ParaView 可视化观察到的问题：

1. 溃坝流体前沿底部区域沿 `x` 方向速度梯度很大。
2. 前沿粒子在 `x` 方向快速拉开，局部最大粒子间距在 `t ≈ 0.2 s` 时接近 `3h`。
3. 最前沿近自由面粒子邻域内严重缺少 `x` 方向邻居。
4. 邻居粒子主要集中在 `y`、`z` 方向，局部邻域呈强各向异性。
5. LSMPS 局部矩阵虽然可能仍有一定邻居数，但几何分布不足以支撑完整三维多项式重构。
6. 算子离散失效后出现异常压力梯度和极大速度，最终导致计算崩溃。

该问题不是单纯的邻居数量不足，而是 **邻域几何质量退化**。即使 `neighbor_count` 没有低到零，只要邻居集中在某些方向，LSMPS moment matrix 仍会病态或失去有效三维重构能力。

## 1.1 阶段一诊断重跑记录

重跑时间：2026-06-23

使用配置：

- 基于 `build/config_dambreak.ini`
- `diagnostics_enabled = false`
- 输出目录：`build/output/dam_break_3d_diagnostics_fast`
- 输出间隔：`0.01 s`

本次运行因单步耗时较高，没有继续等待到 `t ≈ 0.2 s`，而是在输出到 `dam_break_3d_fluid_00004.vtk` 后中止。该阶段仍可用于验证诊断字段是否正确输出，并观察早期粒子分布退化趋势。

早期输出统计如下：

| 文件 | 最近流体距离最小值 | 最大流体邻居距离最大值 | 数密度比最小值 | 几何条件数最大值 | 最大速度 |
|---|---:|---:|---:|---:|---:|
| `fluid_00001` | 0.04989 | 0.15047 | 0.29617 | 2.15994 | 0.23031 |
| `fluid_00002` | 0.04943 | 0.15229 | 0.29250 | 2.20927 | 0.47751 |
| `fluid_00003` | 0.04856 | 0.15431 | 0.28586 | 2.33430 | 0.73218 |
| `fluid_00004` | 0.04725 | 0.15495 | 0.27643 | 2.54634 | 0.98869 |

早期趋势：

1. `nearest_fluid_distance` 持续下降，说明局部已经出现粒子过密趋势。
2. `max_fluid_neighbor_distance` 持续上升，说明局部也在出现拉伸和稀疏趋势。
3. `number_density_ratio` 最小值持续下降，说明最不利区域的有效数密度正在恶化。
4. `geometry_condition_number` 最大值持续上升，说明邻域方向质量正在变差。
5. 最大速度快速增长，说明前沿加速会继续放大粒子分布退化。

当前早期诊断尚未直接达到用户观察到的 `t ≈ 0.2 s` 崩溃阶段，但趋势与可视化观察一致：粒子分布质量正在同时向过密和过稀两个方向恶化，邻域几何条件数也开始上升。

补充说明：

- 尝试使用 `diagnostics_enabled = true` 跑完整 LSMPS 特征值诊断时，计算速度明显变慢，短时间内只能输出到很早期文件。
- 若后续需要自动捕捉 `t ≈ 0.2 s` 附近的失稳，应优先开发失败步自动输出和更低成本的运行期统计，而不是依赖人工等待长时间 VTK 序列。

## 2. 失稳机理分析

### 2.1 前沿拉伸导致 x 方向邻居缺失

溃坝前沿底部速度最大，且沿流向速度差显著。相邻粒子满足：

```text
|u_i - u_j| * dt / h
```

过大时，一个时间步内粒子间距会被明显拉开。连续多个时间步后，前沿粒子在 `x` 方向形成空洞。

当 `x` 方向缺少邻居时，LSMPS 对 `∂/∂x`、`∂²/∂x²` 以及混合空间分布的识别能力下降。局部基函数矩阵不再代表完整三维邻域。

### 2.2 底壁截断放大前沿问题

前沿底部粒子同时受到两个不利因素影响：

1. 自由面或近自由面导致上游/空气侧邻域缺失。
2. 底壁单层壁面导致下方邻域被壁面截断。

因此前沿底部粒子的有效流体邻居本来就少。一旦 `x` 方向再被拉开，邻域会退化成近似二维甚至一维分布。

### 2.3 过稀和过密会同时出现

溃坝前沿不是只出现稀疏区。局部速度差也可能造成粒子追赶和堆积，从而产生过密区。

过稀会导致：

- LSMPS 矩阵秩亏或病态。
- 压力梯度离散失效。
- PPE 局部行质量下降。

过密会导致：

- 权重函数局部过大。
- 局部压力或速度修正异常。
- 粒子穿插后产生更强的非物理相对速度。

因此后续稳定化不能只处理稀疏问题，也必须加入最小粒子间距保护。

## 3. 需要新增的诊断量

当前 VTK 只有总邻居数，无法直接判断方向性退化。建议新增一组稳定性诊断输出。

实现状态：阶段一已完成。当前 TimeStepper 会在流体 VTK 输出中增加粒子分布稳定性诊断字段，用于定位前沿底部粒子在 `x` 方向邻居缺失和邻域几何退化问题。

### 3.1 粒子间距诊断

对每个流体粒子输出：

```text
nearest_fluid_distance
max_fluid_neighbor_distance
mean_fluid_neighbor_distance
number_density
number_density_ratio
```

当前已实现字段：

- `nearest_fluid_distance`
- `max_fluid_neighbor_distance`
- `mean_fluid_neighbor_distance`
- `number_density`
- `number_density_ratio`

用途：

- 判断过密：`nearest_fluid_distance < 0.5h ~ 0.7h`
- 判断过稀：`number_density_ratio` 过低或最大邻居距离过大

### 3.2 方向覆盖诊断

仅有邻居数不足以判断 LSMPS 可用性。建议增加方向覆盖指标。

可选指标：

```text
x_positive_neighbor_count
x_negative_neighbor_count
y_positive_neighbor_count
y_negative_neighbor_count
z_positive_neighbor_count
z_negative_neighbor_count
directional_coverage_score
```

其中方向覆盖分数可以定义为：

```text
coverage = min(count(+x), count(-x), count(+y), count(-y), count(+z), count(-z))
```

或者使用邻居方向张量：

```text
G_i = sum_j w_ij * e_ij * e_ij^T
```

输出：

```text
geometry_min_eigenvalue
geometry_condition_number
```

当前已实现字段：

- `x_positive_neighbor_count`
- `x_negative_neighbor_count`
- `y_positive_neighbor_count`
- `y_negative_neighbor_count`
- `z_positive_neighbor_count`
- `z_negative_neighbor_count`
- `directional_coverage_score`
- `geometry_min_eigenvalue`
- `geometry_max_eigenvalue`
- `geometry_condition_number`

实现位置：

- `src/diagnostics/particle_distribution_diagnostics.hpp`
- `src/diagnostics/particle_distribution_diagnostics.cpp`
- `src/time_integration/time_stepper.cpp`

对应测试：

- `tests/particle_distribution_diagnostics_test.cpp`

该指标能直接识别“邻居都在 y、z 方向，x 方向缺失”的情况。

### 3.3 LSMPS 失败原因诊断

建议细化 LSMPS 状态：

```text
NotEnoughNeighbors
RankDeficient
IllConditioned
DirectionalCoverageFailed
TooSparse
TooDense
```

并输出：

```text
lsmps_regular_status
lsmps_pressure_neumann_status
lsmps_geometry_quality
```

正式长时间模拟可以关闭完整特征值诊断，但稳定性研究阶段应开启完整诊断。

## 4. 时间步控制改进

当前程序的时间步控制已恢复为最原始的 CFL 数约束。也就是只根据邻域内最大相对速度、`max_dt`、`min_dt`、输出时间和结束时间选择时间步，不再使用单步绝对位移限制，也不再在 Correction 阶段对压力修正速度做限幅。

撤回原因：

1. 溃坝模拟在加入单步绝对位移限制和压力修正速度限幅后，仍然在 `t ≈ 0.2 s` 左右崩溃。
2. 对 `build/output/dam_break_3d` 的 VTK 序列统计表明，失稳前主要问题是前沿底部粒子云已经发生明显的方向性退化和局部过密/过稀并存，而不是单一的单步位移过大。
3. 以 `h = 0.05 m`、`displacement_cfl_number = 0.1` 为例，在 `fluid_00021` 附近最大速度约 `4.19 m/s`，绝对位移限制给出的时间步上限约 `0.00119 s`，仍大于当前 `max_dt = 0.001 s`，因此该限制基本没有介入。
4. 以 `pressure_correction_cfl_number = 0.2`、`dt = 0.001 s` 为例，压力修正速度限幅阈值约 `10 m/s`，当前压力修正速度没有触发限幅；该保护只能抑制已经异常放大的速度修正，不能恢复坏邻域几何。

当前代码状态：

- `TimeStepController::chooseDt()` 只接受 `current_time` 和 `max_relative_velocity`。
- `[time]` 配置中不再包含 `displacement_cfl_number` 和 `pressure_correction_cfl_number`。
- `CorrectionParticleStatus` 不再包含 `VelocityCorrectionLimited`。
- Correction 模块直接使用压力梯度得到速度修正，不做额外限幅。

### 4.1 相对速度限制

增加 pairwise relative velocity 限制：

```text
dt_pair = C_pair * h / max_ij(|u_i - u_j|)
```

建议：

```text
C_pair = 0.05 ~ 0.10
```

该限制比单纯最大速度 CFL 更能控制粒子间距拉伸。

### 4.2 单步位移限制

增加最大位移限制：

```text
dt_disp = C_disp * h / max_i(|u_i|)
```

建议：

```text
C_disp = 0.05 ~ 0.10
```

目标是保证单步内粒子不会跨越过大的空间距离。

状态：已试验并撤回。该约束在当前溃坝配置中没有明显改善 `t ≈ 0.2 s` 附近的崩溃问题，原因是粒子云质量退化已经由连续拉伸、局部过密和边界截断共同造成，单纯减小绝对位移不能把缺失的邻域方向补回来。

```text
dt_disp = displacement_cfl_number * h / max_i(|u_i|)
```

最终时间步取：

```text
dt = min(dt_current/growth, dt_cfl, dt_disp, max_dt, output_time, end_time)
```

其中 `dt_cfl` 仍然使用相邻粒子最大相对速度控制：

```text
dt_cfl = cfl_number * h / max_ij(|u_i - u_j|)
```

### 4.3 压力修正速度限制

Correction 后的压力速度修正应限幅：

```text
|delta_u_pressure| <= C_pressure * h / dt
```

这不是物理模型，而是防止局部坏矩阵产生的异常压力梯度把算例瞬间炸掉。

状态：已试验并撤回。该保护属于故障后的被动限幅，不能解决 LSMPS/PPE 之前的邻域几何退化问题。当前阶段优先保留原始压力修正路径，以便后续稳定性改进的效果更清晰。

```text
|delta_u_pressure| <= pressure_correction_cfl_number * h / dt
```

## 4.4 已撤回阶段二试验记录

验证时间：2026-06-23

验证配置：

- 基于 `build/config_dambreak.ini`
- 临时输出目录：`build/output/dam_break_3d_stage2_check`
- `end_time = 0.02`
- `output_interval = 0.01`

输出文件：

- `dam_break_3d_fluid_00000.vtk`
- `dam_break_3d_fluid_00001.vtk`
- `dam_break_3d_fluid_00002.vtk`
- `dam_break_3d_wall_00000.vtk`

短程运行完成 20 个时间步，PPE 均收敛，最大速度从 `0.019949` 增长到 `0.477515`。流体 VTK 已确认包含：

- `nearest_fluid_distance`
- `geometry_condition_number`
- `correction_status`
- `pressure_gradient`

在 `t = 0.02 s` 的输出中，`correction_status` 统计为：

```text
Updated = 3575
VelocityCorrectionLimited = 0
```

后续继续运行表明，加入这两个保护后溃坝模拟仍然在 `t ≈ 0.2 s` 左右崩溃。因此阶段二保护已从当前代码中撤回，后续稳定性改进应转向受约束 particle shifting 和坏粒子保护。

## 5. Particle Shifting 改进方向

普通数密度梯度 shifting 无法直接用于边界附近，因为自由面和壁面本身会造成数密度截断。建议开发 **受约束 particle shifting**。

### 5.1 Shifting 组成

建议第一版使用：

```text
delta_x = delta_x_density + delta_x_repulsion
```

其中：

1. `delta_x_density` 用于缓解局部稀疏和整体不均匀。
2. `delta_x_repulsion` 用于避免粒子间距过小。

### 5.2 数密度均匀化项

对内部粒子计算：

```text
n_i = sum_j W(r_ij)
```

构造平滑位移：

```text
delta_x_density = C_shift * h * F(number_density_gradient)
```

位移必须限幅：

```text
|delta_x_density| <= alpha_shift * h
```

建议：

```text
C_shift = 0.01 ~ 0.03
alpha_shift = 0.02 ~ 0.05
```

### 5.3 最小距离排斥项

当两个流体粒子距离过小时：

```text
r_ij < r_min
```

增加短程排斥：

```text
delta_x_repulsion += C_rep * (r_min - r_ij) * (x_i - x_j) / r_ij
```

建议：

```text
r_min = 0.6h ~ 0.7h
C_rep = 0.1 ~ 0.5
```

排斥项同样必须进入总位移限幅：

```text
|delta_x| <= alpha_total * h
```

### 5.4 粒子类型处理

不同自由面类型应使用不同 shifting 系数：

| 粒子状态 | 建议处理 |
|---|---|
| Internal | 完整 shifting |
| NearFreeSurface | 弱 shifting，系数约 `0.25C_shift` |
| FreeSurface | 默认不 shifting，或只允许切向 shifting |
| Splash | 不 shifting |
| Wall | 不 shifting |

自由面和 splash 不应直接执行普通数密度 shifting，否则会抹平自由面或把粒子推向空气侧。

### 5.5 壁面约束

靠近壁面时，必须投影 shifting 位移。

当前单层壁面粒子有壁面法向，约定法向指向流体区域。对壁面邻居构造约束：

```text
(x_i + delta_x - x_wall) · n_wall >= d_min
```

实现上可以先做近似投影：

```text
if dot(delta_x, n_wall) < 0:
    delta_x = delta_x - dot(delta_x, n_wall) * n_wall
```

对角落处多个壁面法向依次投影，形成可行锥约束。

### 5.6 自由面约束

自由面粒子的外法向可以由自由面识别中的 open direction 或邻域质心估计。若 `n_fs` 指向空气侧，应禁止：

```text
dot(delta_x, n_fs) > 0
```

即不允许 shifting 把粒子推向自由面外。

第一版可以更保守：

```text
FreeSurface: delta_x = 0
NearFreeSurface: 去掉自由面外法向分量并强烈衰减
```

### 5.7 Shifting 后速度处理

第一版建议 shifting 只改位置，不直接改速度。原因是当前 Correction 已经完成速度更新，直接根据 shifting 位移修正速度可能引入额外非物理动量。

如果后续发现质量或动量误差明显，再研究一致性修正：

```text
u_i_new = interpolate(u, x_i + delta_x_i)
```

但这需要额外插值算子和误差验证，不作为第一版内容。

## 6. LSMPS 坏粒子保护

Particle shifting 无法保证所有前沿粒子立刻恢复良好邻域。程序必须在 LSMPS 失败时继续保持数值有界。

### 6.1 几何质量门槛

在构造或使用 LSMPS 矩阵前检查：

```text
neighbor_count >= min_neighbors
directional_coverage_score >= threshold
geometry_condition_number <= threshold
nearest_fluid_distance >= r_min
```

不满足时标记为坏粒子：

```text
BadLsmpsParticle
```

### 6.2 PPE 中的保护

坏粒子不应使用正常压力拉普拉斯离散。可选保护：

1. 将该粒子压力设为局部邻居平均值。
2. 对该粒子使用 Dirichlet-like 保护行。
3. 从 PPE 活跃自由度中临时剔除。

第一版建议采用简单保护行：

```text
p_i = average(p_j)
```

或：

```text
p_i = p_previous_i
```

目标是避免坏粒子生成异常压力尖峰。

### 6.3 Correction 中的保护

对坏粒子的压力梯度进行限幅或置零：

```text
if bad_lsmps_particle:
    grad_p_i = clamp_or_zero(grad_p_i)
```

并限制压力修正速度：

```text
|delta_u_pressure| <= C_u * h / dt
```

## 7. PPE 求解失败保护

当 PPE 不收敛或残差异常时，当前程序直接失败是合理的开发阶段行为。但长时间溃坝模拟需要更细的诊断和保护。

建议输出：

```text
ppe_rhs
ppe_divergence
ppe_wall_neumann_source
ppe_row_status
pressure_limited
velocity_limited
```

正式输出可以保持精简，但失败步应自动输出完整调试 VTK。

## 8. 稳定性开发顺序建议

建议按以下顺序推进。

### 阶段 1：诊断增强

目标：先准确定位何时、何处、因为什么失稳。

状态：已完成第一版实现

任务：

1. 输出最近流体粒子距离。
2. 输出数密度和数密度比。
3. 输出方向覆盖指标。
4. 输出 LSMPS 几何质量指标。
5. 在失败步自动输出调试 VTK。

验收：

- 能在 `t ≈ 0.2 s` 前明确看到前沿底部 `x` 方向覆盖退化。
- 能区分过稀、过密和方向性退化。

### 阶段 2：时间步限制增强

目标：降低粒子分布在单步内被拉坏的概率。

状态：已完成第一版实现

任务：

1. 增加 pairwise 相对速度限制。
2. 增加单步最大位移限制。
3. 增加压力修正速度限幅。

验收：

- 前沿最大粒子间距增长速度降低。
- `t = 0.2 s` 前邻域崩溃明显延后或消失。

### 阶段 3：受约束 Particle Shifting

目标：缓解过稀和过密，同时不破坏壁面和自由面。

任务：

1. 实现数密度均匀化 shifting。
2. 实现最小距离排斥项。
3. 实现壁面法向投影。
4. 实现自由面状态分级处理。
5. 输出 shifting 位移和限幅诊断。

验收：

- 前沿粒子最大间距不再接近 `3h`。
- 最近粒子距离不低于设定阈值。
- 壁面附近粒子不被推入壁面外。
- 自由面不被明显抹平。

### 阶段 4：坏粒子保护

目标：即使局部粒子质量短时间恶化，程序也不会因单个坏粒子崩溃。

任务：

1. 增加 LSMPS 坏粒子状态。
2. PPE 中对坏粒子使用保护行。
3. Correction 中对坏粒子的压力梯度和速度修正限幅。

验收：

- 局部坏粒子不会产生极大压力或速度。
- 溃坝算例能够越过当前 `t ≈ 0.2 s` 的崩溃点。

### 阶段 5：粒子插入或重采样研究

如果前沿区域仍出现不可恢复的稀疏区，则 shifting 不能凭空生成粒子，需要进一步研究粒子插入。

触发条件：

```text
number_density_ratio < threshold
directional_coverage_score too low
nearest gap > 2h
```

可能方案：

1. 在空洞区域插入新粒子。
2. 用邻域插值初始化速度、压力和自由面状态。
3. 对插入粒子设置短期稳定化标记。

该阶段复杂度高，建议在前四个阶段完成后再开始。

## 9. 当前结论

溃坝前沿崩溃的核心原因是流向速度梯度导致粒子在 `x` 方向快速拉开，并与底壁和自由面邻域截断叠加，使前沿近自由面粒子的邻域几何严重各向异性。该问题不能只靠 PPE 或 LSMPS 参数修正解决。

推荐路线是：

1. 先增强诊断，量化粒子间距、方向覆盖和 LSMPS 几何质量。
2. 增强时间步控制，限制相邻粒子相对位移。
3. 开发受约束 particle shifting，同时处理过稀和过密。
4. 增加 LSMPS/PPE/Correction 坏粒子保护，避免局部坏粒子触发全局崩溃。
5. 若仍出现不可恢复稀疏区，再研究粒子插入或重采样。

第一优先级不是直接调大黏性或改 PPE 求解器，而是建立粒子分布质量控制链路。只有邻域几何质量保持可用，LSMPS 离散和 PPE 求解才有稳定基础。
