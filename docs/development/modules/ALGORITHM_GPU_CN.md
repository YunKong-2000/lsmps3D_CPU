# 基于 Cubed Sphere 的三维自由面粒子判定算法与 GPU 实现说明

本文档对应当前实现文件 `surface_detector.cpp`，说明三维粒子法中自由面粒子的几何判定流程，以及后续迁移到 GPU 时建议采用的数据结构和 kernel 设计。

## 1. 问题目标

给定三类信息：

- 流体粒子位置 `x_i`
- 单层壁面粒子位置 `x_w`
- 壁面粒子法向 `n_w`

需要判断每个流体粒子是否为自由面粒子。

这里的壁面法向约定为：

```text
n_w 指向流体侧
```

因此，对某个壁面粒子而言，法向反方向为固体侧。

算法的核心思想是：对每个流体粒子构造一个虚拟球面屏幕，把球面方向分为三类：

```text
solid         该方向进入固体或壁面背后
fluid_shadow 该方向被邻近流体粒子遮挡
open          该方向既不是 solid，也没有被 fluid_shadow 遮挡
```

自由面粒子应该存在较大的、连续的 `open` 区域。内部粒子即使靠近壁面或角落，在扣除 `solid` 区域后，可达区域通常仍被流体邻居遮挡。

## 2. 主要参数

当前代码中的参数位于 `DetectorParams`：

```cpp
double l0 = 0.01;
double screenRadiusFactor = 2.1;
double wallPatchRadiusFactor = 0.85;
double particleRadiusFactor = 0.5;
double openThreshold = 0.22;
double largestThreshold = 0.14;
double coneAngleDeg = 45.0;
double coneThreshold = 0.62;
double minConeAccessibleRatio = 0.40;
int cubedSphereQ = 8;
```

含义如下：

```text
l0                    初始粒子间距
R_screen = 2.1 l0     默认虚拟屏幕半径，也是邻域搜索半径
r_wall = 0.85 l0      单层壁面粒子的局部平面片半径
r_p = 0.5 l0          流体粒子遮挡半径
Q                     cubed sphere 每个立方体面的采样分辨率
N_dir = 6 Q^2         球面方向采样数
gamma = 45 deg        R_cone 使用的锥半角
f_min = 0.40          锥内最小可达面积比例
```

默认 `Q = 8`，因此：

```text
N_dir = 6 * 8 * 8 = 384
```

若程序中的 MPS 邻域统一使用 `3.1 l0`，可通过命令行设置：

```sh
--radius-factor 3.1
```

邻域半径增大后，流体粒子遮挡会增强，`R_open` 和 `R_cone` 的数值分布会变化，阈值需要重新标定。

当前几何测试中，`3.1 l0` 下推荐起始参数为：

```text
T_open = 0.18
T_cone = 0.62
gamma = 45 deg
f_min = 0.40
```

相比 `2.1 l0`，主要调整是降低 `T_open`，因为更大的邻域半径会让自由面粒子的可见开口比例下降。

## 3. Cubed Sphere 球面采样

算法不使用经纬度网格，而使用 cubed sphere。原因是经纬度网格在极点附近面元小、赤道附近面元大，若实现时处理不严谨，会造成方向偏差。Cubed sphere 更适合 GPU 中按规则索引处理。

### 3.1 采样方向

把单位球看作由立方体 6 个面投影得到。每个面划分为 `Q x Q` 个采样点。

对每个面上的局部坐标：

```text
u = -1 + (ix + 0.5) * 2 / Q
v = -1 + (iy + 0.5) * 2 / Q
```

六个面的未归一化方向分别为：

```text
+x 面: ( 1, u, v)
-x 面: (-1, u, v)
+y 面: ( u, 1, v)
-y 面: ( u,-1, v)
+z 面: ( u, v, 1)
-z 面: ( u, v,-1)
```

然后归一化：

```text
omega_m = d / |d|
```

### 3.2 面积权重

普通 gnomonic cubed sphere 不是严格等面积。当前实现为每个方向预计算面积权重：

```text
w_m = du^2 / |d|^3
du = 2 / Q
```

后续计算面积比例时使用加权面积，而不是简单计数。

### 3.3 邻接关系

为了计算最大连续开口区域 `R_largest`，需要在球面采样点之间建立邻接关系。

当前 CPU 实现采用角距离邻接：

```text
dot(omega_i, omega_j) >= cos(2.75 / Q)
```

满足条件的点互为邻居。

GPU 实现中建议把这个邻接表预处理为定长数组，例如：

```text
dirNeighborOffset[m]
dirNeighborCount[m]
dirNeighbors[]
```

如果后续不做 BFS 连通域，而改用锥形开口指标，则可以不保存这个邻接表。

## 4. 背景网格邻域搜索

当前 CPU 实现使用背景网格 `BackgroundGrid` 分别存储流体粒子和壁面粒子。

网格尺寸取：

```text
h_grid = R_screen = 2.1 l0
```

对粒子位置 `x`，网格索引为：

```text
cell = floor(x / h_grid)
```

对每个待检测流体粒子 `i`，搜索其所在单元及周围单元，收集距离满足：

```text
|x_j - x_i| <= R_screen
```

的邻近流体粒子和壁面粒子。

注意：

```text
流体粒子用于生成 fluid_shadow
壁面粒子用于生成 solid
二者不能混用
```

壁面粒子不应当作为普通流体邻居投射阴影，否则贴壁自由面容易漏判。

## 5. 单层壁面粒子的 solid mask

由于壁面由单层粒子表示，没有真实厚度，因此自由面检测时要把每个壁面粒子解释为一个局部平面片。

对壁面粒子 `w`：

```text
位置: x_w
法向: n_w，指向流体侧
局部平面: dot(x - x_w, n_w) = 0
```

对待检测流体粒子 `i` 的某个球面方向 `omega_m`，发射短射线：

```text
x(t) = x_i + t omega_m
0 < t < R_screen
```

若该射线指向壁面背后，则：

```text
den = dot(omega_m, n_w)
den < 0
```

求射线与壁面局部平面的交点参数：

```text
t = dot(x_w - x_i, n_w) / den
```

若：

```text
0 < t < R_screen
```

再计算交点：

```text
x_hit = x_i + t omega_m
```

并判断交点是否落在壁面粒子的局部平面片范围内：

```text
r_tangent = x_hit - x_w - n_w dot(x_hit - x_w, n_w)
|r_tangent| <= r_wall
```

若成立，则该方向为固体方向：

```text
solid[m] = true
```

多壁面交界处直接取并集：

```text
solid = solid_from_wall_1 OR solid_from_wall_2 OR ...
```

这一步是处理壁面角落误判的关键。

## 6. 流体邻居的阴影投射

对流体邻居 `j`：

```text
r_ij = x_j - x_i
d_ij = |r_ij|
e_ij = r_ij / d_ij
```

将流体粒子近似为半径 `r_p = 0.5 l0` 的遮挡球。其在虚拟屏幕方向上的角半径为：

```text
sin(beta_j) = r_p / d_ij
cos(beta_j) = sqrt(1 - (r_p / d_ij)^2)
```

对每个方向 `omega_m`，若：

```text
dot(omega_m, e_ij) >= cos(beta_j)
```

则方向 `omega_m` 被流体邻居遮挡：

```text
fluid_shadow[m] = true
```

但是只对非 solid 方向投射流体阴影：

```text
if !solid[m] and dot(omega_m, e_ij) >= cos(beta_j):
    fluid_shadow[m] = true
```

## 7. 自由面指标

### 7.1 可达区域

从球面总方向中扣除壁面固体方向：

```text
accessible[m] = !solid[m]
```

可达面积为：

```text
A_accessible = sum(w_m for accessible[m])
```

### 7.2 开口区域

既不是固体方向，也没有被流体邻居遮挡的方向为开口：

```text
open[m] = !solid[m] && !fluid_shadow[m]
```

开口面积为：

```text
A_open = sum(w_m for open[m])
```

定义：

```text
R_open = A_open / A_accessible
```

`R_open` 用于衡量总开口面积。

### 7.3 最大连续开口区域

仅使用 `open[m]` 方向，在 cubed sphere 方向邻接图上做连通分量搜索。

对每个 open 连通分量 `C`，计算：

```text
A_C = sum(w_m for m in C)
```

最大连续开口比例为：

```text
R_largest = max(A_C) / A_accessible
```

`R_largest` 用于过滤零散小缝隙。真正自由面通常对应一个较大的连续开口；内部粒子因排列扰动产生的亮区通常较零散。

### 7.4 锥形开口指标 `R_cone`

`R_largest` 需要连通域搜索，CPU 上容易实现，但 GPU 上不够友好。因此当前 CPU 程序同时实现了一个更适合 GPU 的锥形开口指标 `R_cone`，并默认使用 `R_open + R_cone` 进行判定。

对每个球面方向 `q`，定义一个锥形邻域：

```text
cone(q) = { m | dot(omega_q, omega_m) >= cos(gamma) }
```

其中 `gamma` 是锥半角，当前默认：

```text
gamma = 45 deg
```

对每个锥形邻域计算：

```text
A_cone(q) = sum(w_m for m in cone(q))
A_accessible_cone(q) = sum(w_m for m in cone(q) and !solid[m])
A_open_cone(q) = sum(w_m for m in cone(q) and !solid[m] and !fluid_shadow[m])
```

为避免贴近壁面时某个锥几乎全被 `solid` 扣掉，只剩极少数开口方向而导致 `R_cone` 虚高，当前实现要求：

```text
A_accessible_cone(q) / A_cone(q) >= f_min
```

默认：

```text
f_min = 0.40
```

局部锥形开口比例为：

```text
R_cone(q) = A_open_cone(q) / A_accessible_cone(q)
```

取所有方向中的最大值：

```text
R_cone = max_q R_cone(q)
```

`R_cone` 表示：是否存在某个方向附近有一片较大的局部开口。它不是严格连通面积，但能有效区分真实自由面开口和零散噪声亮斑。

### 7.5 最终判定

当前判定条件为：

```text
is_surface = (R_open >= T_open) && (R_cone >= T_cone)
```

默认：

```text
T_open = 0.22
T_cone = 0.62
```

这些阈值应根据实际粒子分布、壁面粒子间距、扰动程度和邻域半径重新标定。

## 8. CPU 参考流程伪代码

```text
预处理:
    生成 cubed sphere 方向 omega[m]
    计算每个方向的面积权重 w[m]
    建立方向邻接表 dir_neighbors[m]

每个时间步:
    清空流体背景网格和壁面背景网格
    将流体粒子插入流体网格
    将壁面粒子插入壁面网格

    for each fluid particle i:
        fluid_neighbors = search fluid grid within R_screen
        wall_neighbors  = search wall grid within R_screen

        solid[m] = false
        shadow[m] = false
        open[m] = false

        for each direction m:
            for each nearby wall particle w:
                if ray from x_i along omega[m] hits wall patch w:
                    solid[m] = true
                    break

        for each fluid neighbor j:
            compute e_ij and cos(beta_j)
            for each direction m:
                if !solid[m] and dot(omega[m], e_ij) >= cos(beta_j):
                    shadow[m] = true

        A_accessible = 0
        A_open = 0
        for each direction m:
            if !solid[m]:
                A_accessible += w[m]
                if !shadow[m]:
                    open[m] = true
                    A_open += w[m]

        R_open = A_open / A_accessible
        R_largest = largest connected weighted open area / A_accessible
        R_cone = max weighted open ratio inside local cones

        is_surface = R_open >= T_open && R_cone >= T_cone
```

## 9. GPU 实现建议

GPU 实现要避免直接照搬 CPU 版本中的动态数组、`std::vector`、`unordered_map` 和 BFS。推荐将流程改为结构化数组和固定长度 mask。

### 9.1 数据布局

建议使用 SoA，而不是 AoS：

```text
fluid_x[Nf], fluid_y[Nf], fluid_z[Nf]
wall_x[Nw], wall_y[Nw], wall_z[Nw]
wall_nx[Nw], wall_ny[Nw], wall_nz[Nw]
```

输出：

```text
is_surface[Nf]
R_open[Nf]
R_largest_or_cone[Nf]
```

Cubed sphere 方向表放入只读缓存或常量内存：

```text
dir_x[Ndir], dir_y[Ndir], dir_z[Ndir]
dir_weight[Ndir]
```

这些方向和权重每次仿真不变，只需初始化一次。

### 9.2 GPU 背景网格

推荐使用排序网格，而不是哈希表链表。

步骤：

```text
1. 对每个流体粒子计算 cell_id
2. 对每个壁面粒子计算 cell_id
3. 分别生成 pair(cell_id, particle_id)
4. 按 cell_id 排序
5. 建立 cell_start[cell] 和 cell_end[cell]
```

可以分别建立：

```text
fluid_cell_start / fluid_cell_end
wall_cell_start  / wall_cell_end
```

如果计算域规则且尺寸已知，`cell_id` 可用三维索引线性化：

```text
cell_id = ix + nx_cell * (iy + ny_cell * iz)
```

如果计算域稀疏或尺寸变化较大，可以使用哈希 cell id，但仍建议排序后用区间索引。

### 9.3 线程映射

最简单、最稳妥的版本：

```text
一个 GPU 线程处理一个流体粒子
```

该线程完成：

```text
邻域遍历
solid mask 计算
fluid shadow 计算
R_open 计算
surface 判定
```

这种实现简单，但每个线程工作量较大。适合先验证。

若粒子数很大或 `Ndir` 较大，可以改为：

```text
一个 warp 或一个 block 处理一个流体粒子
```

方向维度由 warp lanes 或 block 内线程并行处理。

### 9.4 mask 存储

不要在 GPU 线程中使用：

```cpp
bool solid[384];
bool shadow[384];
bool open[384];
```

这种写法会造成大量本地内存访问或寄存器压力。

推荐使用 bit mask：

```text
Ndir = 384
Nword = 384 / 32 = 12

uint32_t solidMask[12]
uint32_t shadowMask[12]
```

方向 `m` 对应：

```text
word = m >> 5
bit  = m & 31
```

置位：

```text
mask[word] |= (1u << bit)
```

开口 mask：

```text
openMask[word] = ~(solidMask[word] | shadowMask[word]) & validMask[word]
```

其中 `validMask` 用于处理 `Ndir` 不是 32 整数倍时的尾部无效 bit。对于 `Q = 8`，`Ndir = 384`，刚好是 32 的整数倍。

### 9.5 `R_open` 的 GPU 计算

若使用 bit mask，计算 `R_open` 有两种方式。

第一种是按方向循环，支持权重：

```text
A_accessible = 0
A_open = 0
for m in 0..Ndir-1:
    if bit(solidMask, m) == 0:
        A_accessible += w[m]
        if bit(shadowMask, m) == 0:
            A_open += w[m]
R_open = A_open / A_accessible
```

第二种是使用 popcount，只适用于近似等面积或不使用权重：

```text
accessible_count += popcount(~solidMask[word])
open_count += popcount(~(solidMask[word] | shadowMask[word]))
R_open = open_count / accessible_count
```

当前 cubed sphere 使用面积权重，因此建议第一版保留权重求和。若后续改用严格等面积 cubed sphere 或 HEALPix，可使用 popcount 加速。

### 9.6 `R_largest` 在 GPU 上的处理

完整 BFS 连通域不适合 GPU 的单线程粒子 kernel，原因是：

```text
分支多
循环次数不固定
需要队列
邻接访问不连续
```

有三种实现路线。

#### 路线 A：先不做 `R_largest`

第一版 GPU 可只计算：

```text
R_open
```

并搭配更严格的阈值或壁面修正密度。优点是实现简单、速度高。

#### 路线 B：用 `R_cone` 替代 `R_largest`

推荐 GPU 版本采用这一方式。

预处理每个方向附近一个固定角度锥内的方向集合：

```text
coneMask[q][word]
coneTotalWeight[q]
```

例如锥角：

```text
gamma = 35 deg ~ 50 deg
```

对每个中心方向 `q`：

```text
open_in_cone = openMask & coneMask[q]
accessible_in_cone = accessibleMask & coneMask[q]
```

计算：

```text
只在 area(accessible_in_cone) / coneTotalWeight[q] >= f_min 时使用该锥
R_cone(q) = area(open_in_cone) / area(accessible_in_cone)
R_cone = max_q R_cone(q)
```

最终判定：

```text
is_surface = R_open >= T_open && R_cone >= T_cone
```

这个指标能近似表达“存在一个连续大开口”，但没有 BFS 那么复杂，更适合 GPU。

#### 路线 C：block 内并行连通域

如果必须保留 `R_largest`，建议：

```text
一个 block 处理一个粒子
```

把 `openMask` 放入 shared memory，在 block 内做固定轮数的连通区域传播或并查集。该方法实现复杂，且通常不如 `R_cone` 性价比高。

### 9.7 推荐 GPU kernel 结构

推荐拆成多个 kernel：

```text
Kernel 1: computeFluidCellId
Kernel 2: sort fluid particles by cell id
Kernel 3: build fluid cell start/end

Kernel 4: computeWallCellId
Kernel 5: sort wall particles by cell id
Kernel 6: build wall cell start/end

Kernel 7: classifySurfaceParticles
```

其中 `classifySurfaceParticles`：

```text
thread i -> fluid particle i

1. 初始化 solidMask/shadowMask = 0
2. 遍历附近 wall cell
   对候选壁面粒子和方向，设置 solidMask
3. 遍历附近 fluid cell
   对候选流体粒子和方向，设置 shadowMask
4. 计算 R_open
5. 计算 R_cone 或简化的连续开口指标
6. 写出 is_surface[i], R_open[i], R_cone[i]
```

### 9.8 性能注意事项

1. 不要在 kernel 内动态分配内存。
2. 不要在每个线程内重新生成 cubed sphere 方向。
3. 不要使用 GPU 上的递归、动态队列或不定长 `vector`。
4. 邻域粒子数量可能差异很大，容易造成 warp divergence。
5. 对 `Q = 8`，`Ndir = 384`，每粒子每邻居都遍历全部方向会偏重。可以预处理方向加速结构。

可选优化：

```text
1. 对每个邻居方向 e_ij，只检查其附近角帽内的方向，而不是全部 Ndir。
2. 为每个方向预存 coneMask，用 bit 操作计算 R_cone。
3. 将 wall solid 判断改为先筛选可能命中的壁面方向范围。
4. 将常用参数放入 constant memory。
5. 使用 half 或 float 存储方向与权重，但几何判断建议先用 float 验证误差。
```

### 9.9 适合 GPU 的简化版判定

当前 CPU 版本已经实现：

```text
Cubed sphere + solidMask + shadowMask + R_open + R_cone
```

并保留 `R_largest` 作为诊断输出。工程上建议 GPU 第一版也使用：

```text
Cubed sphere + solidMask + shadowMask + R_open + R_cone
```

而不是：

```text
Cubed sphere + solidMask + shadowMask + R_open + BFS R_largest
```

推荐判定：

```text
is_surface = R_open >= T_open && R_cone >= T_cone
```

这样可以保留壁面修正和连续开口判定的主要效果，同时避免 GPU 上低效的连通域搜索。

## 10. VTK 可视化输出

当前程序会把流体粒子和壁面粒子同时写入 VTK。

字段：

```text
particle_type        0 为流体，1 为壁面
is_surface           判定结果
expected_surface     几何测试算例中的期望标签
has_expected         是否有期望标签
R_open               总开口比例
R_largest            最大连续开口比例
classification_code  分类结果编码
wall_normal          壁面法向，流体粒子为零向量
```

`classification_code`：

```text
0 未标注
1 TP，正确判为自由面
2 TN，正确判为内部
3 FP，内部误判为自由面
4 FN，自由面漏判为内部
```

在 ParaView 中建议：

```text
1. 用 particle_type 区分流体和壁面
2. 用 is_surface 查看自由面判定
3. 用 classification_code 查看误判和漏判
4. 用 R_open/R_largest 检查阈值是否合理
```

## 11. 当前测试算例

当前实现提供四个内置算例：

```text
hydrostatic
dam-break
large-hydrostatic
large-dam-break
```

其中大算例可用：

```sh
./surface_detector --case large-hydrostatic --l0 0.02 --q 8 \
  --output large_hydrostatic.csv --vtk large_hydrostatic.vtk

./surface_detector --case large-dam-break --l0 0.02 --q 8 \
  --output large_dam_break.csv --vtk large_dam_break.vtk
```

这两个算例用于几何判定验证，不包含完整流体动力学时间推进。若要验证真实 MPS 计算中的鲁棒性，应读取实际时刻的粒子云，并用人工标注、界面重建结果或高分辨率参考结果作为期望标签。
