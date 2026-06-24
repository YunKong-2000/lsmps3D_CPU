# Particle Shifting 最小实现方案

本文档记录当前程序中 Particle Shifting，简称 PS，模块的最小实现方案。当前阶段只实现**流体粒子之间的短程排斥项**，用于避免局部粒子间距过小。

当前版本不实现：

1. 拉伸回拉项。
2. 数密度梯度项。
3. 壁面投影。
4. 自由面投影。
5. 速度重映射。
6. 粒子插入或重采样。

自由面粒子和 splash 粒子暂时不做 PS。

## 1. 当前目标

溃坝模拟中前沿区域同时存在两类粒子分布问题：

1. 流向拉伸导致粒子间距过大。
2. 局部追赶或穿插导致粒子间距过小。

当前最小实现只处理第二类问题，即当两个可移动流体粒子距离小于给定阈值时，对粒子施加小幅排斥位移，避免粒子继续靠得过近。

## 2. 配置参数

当前只开放三个参数：

```ini
[particle_shifting]
enabled = true
max_displacement_factor = 0.05
min_distance_factor = 0.70
```

对应数学量：

$$
\alpha = \text{max\_displacement\_factor}
$$

$$
r_{min} = \text{min\_distance\_factor} \cdot h
$$

其中 $h$ 是初始粒子间距。

默认值为：

$$
\alpha = 0.05
$$

$$
r_{min}=0.70h
$$

因此对 $h=0.05m$ 的溃坝算例，单步最大 PS 位移为：

$$
\alpha h = 0.0025m
$$

## 3. 作用对象

当前只移动以下流体粒子：

| 粒子状态 | 是否执行 PS |
|---|---|
| `Internal` | 是 |
| `NearFreeSurface` | 是 |
| `FreeSurface` | 否 |
| `Splash` | 否 |
| `Wall` | 否 |

其中 `FreeSurface` 和 `Splash` 暂时完全跳过，避免第一版 PS 直接改变自由面形态。

## 4. 排斥项公式

对流体粒子 $i$，遍历其流体邻居 $j \in \mathcal{N}_i^f$。两粒子距离为：

$$
r_{ij}
=
\left\|
\mathbf{x}_i-\mathbf{x}_j
\right\|
$$

若：

$$
r_{ij} \ge r_{min}
$$

则该邻居不产生 PS 贡献。

若：

$$
r_{ij} < r_{min}
$$

则计算排斥方向：

$$
\mathbf{d}_{ij}
=
\frac{
\mathbf{x}_i-\mathbf{x}_j
}{
r_{ij}+\epsilon
}
$$

排斥强度：

$$
s_{ij}
=
\frac{r_{min}-r_{ij}}{r_{min}}
$$

粒子 $i$ 的未限幅排斥方向和为：

$$
\mathbf{q}_i
=
\sum_{j \in \mathcal{N}_i^f,\ r_{ij}<r_{min}}
s_{ij}\mathbf{d}_{ij}
$$

原始 PS 位移为：

$$
\Delta \mathbf{x}_i^{raw}
=
\alpha h \mathbf{q}_i
$$

最后进行单步位移限幅：

$$
\Delta \mathbf{x}_i^{ps}
=
\operatorname{clip}
\left(
\Delta \mathbf{x}_i^{raw},
\alpha h
\right)
$$

其中：

$$
\left\|
\Delta \mathbf{x}_i^{ps}
\right\|
\le
\alpha h
$$

## 5. 时间步接入位置

当前 PS 放在 Correction 之后：

1. 邻域搜索。
2. 自由面识别。
3. LSMPS 矩阵构造。
4. Provisional 速度计算。
5. PPE 求解。
6. Correction 更新速度和位置。
7. Particle Shifting 对位置做短程排斥修正。
8. 下一时间步重新邻域搜索。

PS 只修改位置，不修改速度：

$$
\mathbf{x}_i^{n+1}
\leftarrow
\mathbf{x}_i^{n+1}
+
\Delta \mathbf{x}_i^{ps}
$$

$$
\mathbf{u}_i^{n+1}
\leftarrow
\mathbf{u}_i^{n+1}
$$

当前这么处理是为了把 PS 定义为点集质量修正，而不是额外的物理速度。

## 6. VTK 诊断输出

当 `[particle_shifting] enabled = true` 时，流体 VTK 增加：

向量场：

```text
particle_shift
```

标量场：

```text
particle_shift_magnitude
particle_shift_limited
particle_shift_repulsion_active
```

这些字段用于检查：

1. 哪些粒子触发了排斥。
2. PS 位移大小是否接近上限。
3. 排斥是否主要发生在过密区域。

## 7. 当前实现位置

代码位置：

- `src/particle_shifting/particle_shifter.hpp`
- `src/particle_shifting/particle_shifter.cpp`

配置：

- `SimulationConfig::particle_shifting`
- `[particle_shifting]` INI 配置段

测试：

- `tests/particle_shifter_test.cpp`

TimeStepper 接入：

- Correction 后计算并应用 PS。
- 输出 VTK 时写出 PS 诊断字段。

## 8. 当前限制

当前最小实现只能防止粒子过近，不能解决前沿粒子被拉开导致的空洞问题。因此它可能缓解局部过密导致的不稳定，但不能单独解决溃坝前沿接近 $3h$ 的粒子间距拉伸。

如果后续验证表明前沿拉伸仍然是主要崩溃原因，再逐步加入拉伸回拉项或其他粒子分布修正方法。
