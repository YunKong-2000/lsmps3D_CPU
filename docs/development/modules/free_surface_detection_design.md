# 自由面粒子识别模块实现说明

本文档记录当前自由面识别模块的实际实现路径。当前实现参考 `docs/development/modules/ALGORITHM_GPU_CN.md` 中的 GPU 友好算法，采用 cubed sphere 球面采样、单层壁面粒子的 solid mask、流体邻居的 shadow mask，并使用 `R_open + R_cone` 综合判定自由面。

当前版本不实现最大连通域判定 `R_largest`。

## 1. 代码位置

核心实现：

- `src/free_surface/free_surface_detector.hpp`
- `src/free_surface/free_surface_detector.cpp`

相关依赖：

- `src/core/particle_set.*`
- `src/core/particle_types.hpp`
- `src/core/simulation_config.*`
- `src/neighbor/neighbor_search.*`
- `src/io/vtk_writer.*`

测试与可视化验收：

- `tests/free_surface_detector_test.cpp`
- `output/free_surface_open_top.vtk`
- `output/free_surface_closed_box.vtk`
- `output/free_surface_wall_contact.vtk`
- `output/free_surface_splash.vtk`

## 2. 状态枚举

流体粒子状态定义在 `src/core/particle_types.hpp`：

```cpp
enum class FluidParticleState {
    Internal = 0,
    NearFreeSurface = 1,
    FreeSurface = 2,
    Splash = 3,
};
```

VTK 输出中的 `fluid_state` 使用同样的整数映射：

```text
0 = Internal
1 = NearFreeSurface
2 = FreeSurface
3 = Splash
```

## 3. 输入与输出

### 3.1 输入

`FreeSurfaceDetector::detect()` 的输入为：

```cpp
FreeSurfaceDiagnostics detect(
    ParticleSet& particles,
    const TypedNeighborList& neighbors) const;
```

其中：

- `particles.positions()` 提供粒子坐标。
- `particles.types()` 区分流体粒子和壁面粒子。
- `particles.wallNormals()` 提供壁面粒子法向，约定为指向流体侧。
- `neighbors.fluid[i]` 是粒子 `i` 的流体邻居索引。
- `neighbors.wall[i]` 是粒子 `i` 的壁面邻居索引。

如果壁面粒子法向为零向量，当前实现会跳过该壁面粒子的 solid mask 贡献。因此后续几何/壁面生成模块必须正确写入壁面法向。

### 3.2 输出

直接写入：

- `particles.fluidStates()`

返回诊断：

```cpp
struct FreeSurfaceDiagnostics {
    std::vector<double> open_ratio;
    std::vector<double> cone_ratio;
    std::vector<double> accessible_area_ratio;
    std::vector<int> reason_code;
};
```

当前 `reason_code` 含义：

```text
0 = Internal 或未触发特殊原因
1 = FreeSurface by R_open + R_cone
3 = NearFreeSurface by distance threshold
4 = Splash
5 = No accessible area
```

`reason_code = 2`、`6` 等在早期设计中预留，但当前实现没有使用。

## 4. 配置参数

配置结构在 `SimulationConfig::free_surface` 中：

```cpp
struct FreeSurfaceConfig {
    double neighbor_count_ratio = 0.85;
    double number_density_ratio = 0.90;
    std::size_t near_surface_layers = 1;
    double screen_radius_factor = 3.1;
    double wall_patch_radius_factor = 0.85;
    double particle_radius_factor = 0.5;
    double open_threshold = 0.18;
    double cone_angle_degrees = 45.0;
    double cone_threshold = 0.62;
    double min_cone_accessible_ratio = 0.40;
    std::size_t cubed_sphere_q = 8;
    std::size_t splash_max_fluid_neighbors = 4;
    double splash_open_threshold = 0.75;
    double near_surface_distance_factor = 1.5;
};
```

示例 INI 配置：

```ini
[free_surface]
screen_radius_factor = 3.1
wall_patch_radius_factor = 0.85
particle_radius_factor = 0.5
open_threshold = 0.18
cone_angle_degrees = 45.0
cone_threshold = 0.62
min_cone_accessible_ratio = 0.40
cubed_sphere_q = 8
splash_max_fluid_neighbors = 4
splash_open_threshold = 0.75
near_surface_distance_factor = 1.5
```

当前 detector 使用的实际半径为：

```text
R_screen = screen_radius_factor * particle_spacing
r_wall   = wall_patch_radius_factor * particle_spacing
r_p      = particle_radius_factor * particle_spacing
```

## 5. 初始化流程

`FreeSurfaceDetector` 构造时完成两个预处理。

### 5.1 Cubed sphere 方向表

使用 `cubed_sphere_q = Q`，生成：

```text
N_dir = 6 * Q * Q
```

默认 `Q = 8`，所以：

```text
N_dir = 384
```

每个方向保存：

- 单位方向 `directions_[m]`
- 面积权重 `direction_weights_[m]`

权重计算为：

```text
du = 2 / Q
w = du^2 / |d|^3
```

其中 `d` 是 cubed sphere 面上的未归一化方向。

### 5.2 Cone 查表

根据 `cone_angle_degrees` 预先为每个方向 `q` 建立锥形邻域：

```text
cone(q) = { m | dot(direction[q], direction[m]) >= cos(cone_angle) }
```

同时预计算每个 cone 的总面积：

```text
cone_area[q] = sum(weight[m] for m in cone(q))
```

这样每次识别时只需遍历查表结果，不需要重复做方向集合构建。

## 6. 单粒子识别流程

当前 CPU 代码每次处理一个流体粒子。这个流程天然可映射为 GPU 中“一个线程处理一个流体粒子”的第一版实现。

对流体粒子 `i`：

```text
solid[m]  = false
shadow[m] = false
```

其中 `m` 是 cubed sphere 方向索引。

### 6.1 壁面 solid mask

对每个壁面邻居 `w`：

```text
x_w = wall position
n_w = wall normal, pointing to fluid side
```

对每个方向 `omega_m` 从流体粒子发射短射线：

```text
x(t) = x_i + t * omega_m
```

若射线朝向壁面背后：

```text
den = dot(omega_m, n_w)
den < 0
```

求与壁面局部平面的交点：

```text
t = dot(x_w - x_i, n_w) / den
```

若：

```text
0 < t < R_screen
```

并且交点落在壁面粒子的局部平面片内：

```text
|tangent_offset| <= r_wall
```

则：

```text
solid[m] = true
```

多个壁面粒子的 solid 区域取并集。

这一部分是当前实现处理单层壁面粒子和壁面角落误判的核心。

### 6.2 流体 shadow mask

对每个流体邻居 `j`：

```text
r_ij = x_j - x_i
d_ij = |r_ij|
e_ij = r_ij / d_ij
```

把流体邻居近似为半径：

```text
r_p = particle_radius_factor * particle_spacing
```

的遮挡球。角半径满足：

```text
cos(beta) = sqrt(1 - min(1, r_p / d_ij)^2)
```

对每个方向 `omega_m`：

```text
if !solid[m] and dot(omega_m, e_ij) >= cos(beta):
    shadow[m] = true
```

注意：当前实现只在非 solid 方向上投射流体阴影，避免壁面背后方向影响自由面开口判断。

## 7. R_open 与 R_cone

### 7.1 可达面积

```text
accessible[m] = !solid[m]
A_accessible = sum(w[m] for accessible[m])
```

如果 `A_accessible` 近似为零，该粒子不参与自由面判定，`reason_code = 5`。

### 7.2 开口面积与 R_open

```text
open[m] = !solid[m] && !shadow[m]
A_open = sum(w[m] for open[m])
R_open = A_open / A_accessible
```

`R_open` 衡量可达球面上没有被流体遮挡的总开口比例。

### 7.3 锥形开口 R_cone

对每个预计算锥 `cone(q)`：

```text
A_cone_accessible = sum(w[m] for m in cone(q) and !solid[m])
A_cone_open       = sum(w[m] for m in cone(q) and !solid[m] and !shadow[m])
```

为避免贴壁时锥内可达面积太小导致比例虚高，先检查：

```text
A_cone_accessible / A_cone >= min_cone_accessible_ratio
```

满足后计算：

```text
R_cone(q) = A_cone_open / A_cone_accessible
```

最终：

```text
R_cone = max_q R_cone(q)
```

当前实现不计算最大连通开口区域 `R_largest`。

## 8. Primary 状态判定

当前第一轮只产生临时状态 `primary_states`：

```text
Internal
FreeSurface
Splash
```

### 8.1 Splash

`Splash` 优先于普通自由面判定：

```text
if fluid_neighbor_count <= splash_max_fluid_neighbors
   and wall_neighbor_count == 0
   and R_open >= splash_open_threshold:
       Splash
```

默认：

```text
splash_max_fluid_neighbors = 4
splash_open_threshold = 0.75
```

该判据刻意要求没有壁面邻居，避免壁面附近稀疏粒子被误判为飞溅粒子。

### 8.2 FreeSurface

若未被判为 `Splash`，再判断普通自由面：

```text
if R_open >= open_threshold and R_cone >= cone_threshold:
    FreeSurface
```

默认：

```text
open_threshold = 0.18
cone_threshold = 0.62
```

## 9. NearFreeSurface 判定

`NearFreeSurface` 在第二轮生成。

对非 `FreeSurface`、非 `Splash` 的流体粒子 `i`，遍历其流体邻居。如果存在 `primary_states[j] == FreeSurface`，并且：

```text
|x_j - x_i| < near_surface_distance_factor * particle_spacing
```

则：

```text
state[i] = NearFreeSurface
reason_code[i] = 3
```

否则：

```text
state[i] = Internal
```

默认：

```text
near_surface_distance_factor = 1.5
```

这里使用严格小于 `<`。也就是说，若阈值为 `1.5h`，距离为 `h` 的粒子会成为近自由面粒子，距离为 `2h` 的粒子不会。

`Splash` 粒子不会诱导周围粒子成为 `NearFreeSurface`。

## 10. VTK 诊断输出

自由面测试会输出以下文件：

```text
output/free_surface_open_top.vtk
output/free_surface_closed_box.vtk
output/free_surface_wall_contact.vtk
output/free_surface_splash.vtk
```

每个文件包含默认粒子字段：

```text
particle_type
fluid_state
neighbor_count
fluid_neighbor_count
wall_neighbor_count
```

以及自由面诊断字段：

```text
free_surface_open_ratio
free_surface_cone_ratio
free_surface_accessible_area_ratio
free_surface_reason_code
```

这些字段用于 ParaView 中检查：

- 自由面是否连续。
- 壁面附近是否误判。
- Splash 是否只出现在孤立小团粒子中。
- `R_open` 和 `R_cone` 的阈值是否合理。

## 11. 当前测试覆盖

测试文件：

```text
tests/free_surface_detector_test.cpp
```

当前覆盖：

1. 开顶流体块：
   - 顶层中心粒子应为 `FreeSurface`。
   - 内部粒子不能为 `FreeSurface`。

2. 封闭盒：
   - 中心内部粒子不能为 `FreeSurface`。
   - 靠壁内部粒子不能为 `FreeSurface`。

3. 靠壁自由面：
   - 自由面与侧壁相交时，靠壁顶部粒子应为 `FreeSurface`。
   - 侧壁中部内部粒子不能为 `FreeSurface`。

4. Splash：
   - 孤立小团粒子中目标粒子应为 `Splash`。
   - `reason_code = 4`。

5. 近自由面距离阈值：
   - 显式检查 `NearFreeSurface = 1`、`FreeSurface = 2`。
   - 将 `near_surface_distance_factor` 降为 `0.5` 时，距离一层粒子间距的粒子不会被标记为 `NearFreeSurface`。

## 12. GPU 迁移路径

当前 CPU 实现仍使用：

```text
std::vector<unsigned char> solid
std::vector<unsigned char> shadow
std::vector<std::vector<std::size_t>> neighbor list
```

后续迁移到 GPU 时建议改为：

```text
fluid_neighbor_offsets[N + 1]
fluid_neighbor_indices[Mf]
wall_neighbor_offsets[N + 1]
wall_neighbor_indices[Mw]
solid_mask bitset
shadow_mask bitset
```

推荐 kernel/pass 划分：

1. 每个流体粒子构造 solid/shadow mask，并计算 `R_open`、`R_cone`。
2. 根据 `R_open`、`R_cone` 和 Splash 判据生成 `primary_state`。
3. 根据最近自由面距离生成 `NearFreeSurface`。
4. 可选加入时间滞回，抑制状态抖动。

当前实现没有使用 BFS 或最大连通域搜索，因此比 `R_largest` 判据更适合 GPU 并行化。

## 13. 当前限制

- 需要壁面粒子具有可靠法向；法向为零的壁面粒子会被跳过。
- 当前没有时间滞回，自由面状态可能在动态算例中抖动。
- 当前没有最大连通域过滤，`R_cone` 用于替代该功能。
- 当前没有直接处理内部局部稀疏区的高级抑制，只依赖 `R_open + R_cone` 和壁面 solid mask。
- 参数仍需要通过静水箱、溃坝和晃荡算例继续标定。
