# 静水算例性能基准与优化节点分析

本文档记录当前 CPU 版本在约 8000 个流体粒子的静水算例中的单步耗时，用于确定后续加速工作的优先级。

## 0. 最新优化成果汇总

最新一次运行日期：2026-06-23

最新运行配置：

- 流体粒子数：7938
- 壁面粒子数：2377
- 总粒子数：10315
- 时间步数：2
- `number_density_ratio = 0.98`
- `lsmps.diagnostics_enabled = false`
- 不采用跨时间步 LSMPS 矩阵复用
- 不采用邻居几何 stencil 缓存
- 暂不采用 OpenMP 并行

最新一次静水算例的原始输出：

```text
[timing] step=1 time=0.001 neighbor=0.163212 time_control=0.0155237 free_surface=1.70161 lsmps=4.51146 provisional=1.29959 ppe=6.94685 correction=0.452294 diagnostics=0.0157623 vtk=0.080858 total=15.1872
[timing] step=2 time=0.002 neighbor=0.161024 time_control=0.0152836 free_surface=1.72168 lsmps=4.58467 provisional=1.31661 ppe=6.99339 correction=0.454852 diagnostics=0.0154167 vtk=0.0741251 total=15.337
```

最新平均耗时统计：

| 模块 | 平均耗时 s | 占总耗时 |
|---|---:|---:|
| PPE 构造与求解 | 6.97012 | 45.67% |
| LSMPS 矩阵构造 | 4.54806 | 29.80% |
| 自由面识别 | 1.71165 | 11.22% |
| Provisional | 1.30810 | 8.57% |
| Correction | 0.45357 | 2.97% |
| 邻域搜索 | 0.16212 | 1.06% |
| VTK 输出 | 0.07749 | 0.51% |
| 时间步控制 | 0.01540 | 0.10% |
| 诊断统计 | 0.01559 | 0.10% |
| 总耗时 | 15.26210 | 100% |

相对最初基准的累计优化效果：

| 指标 | 最初基准 s | 最新结果 s | 下降幅度 |
|---|---:|---:|---:|
| 总单步耗时 | 21.17425 | 15.26210 | 27.92% |
| PPE 构造与求解 | 9.64144 | 6.97012 | 27.71% |
| LSMPS 矩阵构造 | 7.05567 | 4.54806 | 35.54% |
| 自由面识别 | 2.50681 | 1.71165 | 31.72% |

分阶段优化效果：

| 阶段 | 主要优化 | 平均总耗时 s | 相对上一阶段下降 |
|---|---|---:|---:|
| 初始基准 | 未优化 | 21.17425 | - |
| 第一轮 | PPE 散度直接读取邻域临时速度 | 18.64495 | 11.95% |
| 第二轮 | 无壁面邻居粒子跳过 `pressure_neumann` 重复构造 | 17.03455 | 8.64% |
| 第三轮 | PPE 按行批量装配、LSMPS 诊断开关、自由面数密度粗筛 | 15.26210 | 10.41% |

当前瓶颈排序已经变为：

1. PPE 构造与求解，约 45.67%。
2. LSMPS 矩阵构造，约 29.80%。
3. 自由面识别，约 11.22%。
4. Provisional，约 8.57%。

当前结论：已有优化使约 8000 流体粒子的静水算例单步耗时从 `21.17425 s` 降至 `15.26210 s`，累计下降约 `27.92%`。LSMPS 和自由面识别的优化效果明显；PPE 仍是当前最大瓶颈，后续若继续加速，应优先评估 PETSc KSP/PC 参数和求解器路径。

## 1. 基准算例

算例目录：

- `cases/hydrostatic_box`

前处理程序：

- `cases/hydrostatic_box/preprocess_hydrostatic_box.cpp`

运行命令：

```bash
cmake --build build
./build/hydrostatic_box_preprocess
./build/lsmps3d cases/hydrostatic_box/config.ini
```

几何与粒子规模：

- 水箱尺寸：1.0 m x 1.0 m x 1.0 m
- 初始水深：0.85 m
- 粒子间距：0.05 m
- 支持域半径：0.155 m，即 3.1 倍粒子间距
- 流体粒子数：7938
- 壁面粒子数：2377
- 总粒子数：10315

时间设置：

- `initial_dt = 0.001`
- `end_time = 0.002`
- 共运行 2 个时间步
- 每步输出 VTK，并打印模块耗时

## 2. 原始耗时输出

第 1 步：

```text
[timing] step=1 time=0.001 neighbor=0.154496 time_control=0.014587 free_surface=2.49307 lsmps=6.94448 provisional=1.24662 ppe=9.60492 correction=0.443512 diagnostics=0.0154455 vtk=0.067598 total=20.9847
```

第 2 步：

```text
[timing] step=2 time=0.002 neighbor=0.161799 time_control=0.0156626 free_surface=2.52055 lsmps=7.16686 provisional=1.28623 ppe=9.67796 correction=0.444633 diagnostics=0.0156522 vtk=0.0744813 total=21.3638
```

## 3. 平均耗时统计

| 模块 | 平均耗时 s | 占总耗时 |
|---|---:|---:|
| PPE 构造与求解 | 9.64144 | 45.53% |
| LSMPS 矩阵构造 | 7.05567 | 33.32% |
| 自由面识别 | 2.50681 | 11.84% |
| Provisional | 1.26643 | 5.98% |
| Correction | 0.44407 | 2.10% |
| 邻域搜索 | 0.15815 | 0.75% |
| VTK 输出 | 0.07104 | 0.34% |
| 时间步控制 | 0.01512 | 0.07% |
| 诊断统计 | 0.01555 | 0.07% |
| 总耗时 | 21.17425 | 100% |

当前最主要瓶颈排序：

1. PPE 构造与求解
2. LSMPS 矩阵构造
3. 自由面识别
4. Provisional

邻域搜索和 VTK 输出在当前规模下不是主要瓶颈。

## 4. PPE 瓶颈分析

相关代码：

- `src/pressure_poisson/pressure_poisson.cpp`

当前 PPE 平均耗时约 9.64 s，占总时间约 45.5%。主要原因包括：

- 每个时间步都重新装配 PETSc AIJ 矩阵。
- 当前 PPE 每个非自由面粒子都会计算速度散度。
- `computeDivergence` 内部对每个粒子重复构造 `ux`、`uy`、`uz` 三个全粒子数组，这会产生明显的重复开销。
- 当前求解器默认使用 GMRES + ILU，静水算例中约 23 次迭代。
- 目前单进程运行，没有利用 PETSc MPI 并行装配和求解。

优先优化建议：

1. 已完成：优化 `computeDivergence`，避免每个粒子重复构造全场 `ux/uy/uz` 数组。
2. 将速度分量数组在 PPE 外层预构造一次，或直接按邻居索引读取 `provisional.provisional_velocity[j]`。
3. 对静水和小变形算例，考虑复用上一时间步 PPE 矩阵结构，只更新数值和 RHS。
4. 研究更合适的 PETSc 预条件器，例如 BJACOBI+ILU、ASM、HYPRE。
5. 中期实现 PETSc 多进程装配和求解。

### 4.1 已完成优化：PPE 散度计算直接读取邻域速度

优化时间：2026-06-23

优化位置：

- `src/pressure_poisson/pressure_poisson.cpp`

优化前，`computeDivergence` 在每个非自由面流体粒子上都会重新分配并填充三个全粒子数组：

- `ux`
- `uy`
- `uz`

这意味着约 8000 个流体粒子的基准算例中，每个时间步会重复执行约 8000 次全场速度分量数组构造。每个粒子的散度计算实际只需要读取其邻域粒子的速度，因此该全场数组构造属于不必要的重复开销。

优化后，PPE 散度计算不再构造 `ux/uy/uz` 全局分量数组，而是在邻居循环中直接读取：

- `provisional.provisional_velocity[j]`
- `provisional.provisional_velocity[i]`

并分别累加 x、y、z 速度分量对应的 LSMPS 右端项。该实现更接近后续 GPU 并行路径：每个粒子的散度计算只依赖自身邻域数据。

优化后重新运行同一静水基准算例，原始耗时输出为：

第 1 步：

```text
[timing] step=1 time=0.001 neighbor=0.160815 time_control=0.0150163 free_surface=2.50195 lsmps=7.11384 provisional=1.27794 ppe=6.91803 correction=0.44896 diagnostics=0.0157799 vtk=0.0683557 total=18.5207
```

第 2 步：

```text
[timing] step=2 time=0.002 neighbor=0.16569 time_control=0.0154878 free_surface=2.54235 lsmps=7.23444 provisional=1.30099 ppe=6.96788 correction=0.453188 diagnostics=0.0157428 vtk=0.0733816 total=18.7692
```

优化后平均耗时统计：

| 模块 | 平均耗时 s | 占总耗时 |
|---|---:|---:|
| LSMPS 矩阵构造 | 7.17414 | 38.48% |
| PPE 构造与求解 | 6.94295 | 37.24% |
| 自由面识别 | 2.52215 | 13.53% |
| Provisional | 1.28947 | 6.92% |
| Correction | 0.45107 | 2.42% |
| 邻域搜索 | 0.16325 | 0.88% |
| VTK 输出 | 0.07087 | 0.38% |
| 时间步控制 | 0.01525 | 0.08% |
| 诊断统计 | 0.01576 | 0.08% |
| 总耗时 | 18.64495 | 100% |

与优化前相比：

- PPE 平均耗时：`9.64144 s` -> `6.94295 s`，下降约 `27.99%`。
- 总单步耗时：`21.17425 s` -> `18.64495 s`，下降约 `11.95%`。
- 当前最大瓶颈转为 LSMPS 矩阵构造，PPE 仍然是第二大瓶颈。

## 5. LSMPS 瓶颈分析

相关代码：

- `src/lsmps/lsmps_matrices.cpp`

当前 LSMPS 平均耗时约 7.06 s，占总时间约 33.3%。主要原因包括：

- 每个时间步对所有流体粒子重新构造 `regular` 和 `pressure_neumann` 两套逆矩阵。
- 每套矩阵都会做特征值分解、秩判断、条件数估计和 LU 求逆。
- 对无壁面邻居的内部粒子，`pressure_neumann` 与 `regular` 的构造高度重复。
- 为保持 LSMPS 离散精度，当前约定每个时间步必须重新计算 LSMPS 矩阵，不采用跨时间步矩阵复用优化。

优先优化建议：

1. 已完成：对无壁面邻居粒子，直接令 `pressure_neumann = regular`，避免重复构造第二套矩阵。
2. 已完成：增加 `lsmps.diagnostics_enabled` 配置项，允许正式模拟跳过全量特征值诊断和条件数估计。
3. 暂不采用邻居几何 stencil 缓存，避免显著增加内存占用。
4. OpenMP 粒子并行暂不实现，保留为后续方向。

### 5.1 已完成优化：无壁面邻居粒子跳过 pressure_neumann 构造

优化时间：2026-06-23

优化位置：

- `src/lsmps/lsmps_matrices.cpp`

优化前，每个流体粒子都会构造两套 LSMPS 逆矩阵：

- `regular`
- `pressure_neumann`

其中 `pressure_neumann` 只有在邻域中存在壁面粒子时才和 `regular` 不同。对于邻域中没有壁面粒子的内部粒子，两套矩阵的 raw moment matrix 完全相同，因此重复构造 `pressure_neumann` 会额外执行一次邻居遍历、特征值诊断和 LU 求逆。

优化后，如果 `neighbors.wall[i].empty()`，则直接执行：

```cpp
matrices.particles[i].pressure_neumann = matrices.particles[i].regular;
```

该优化不改变 LSMPS 离散结果，也不改变后续 PPE 和 Correction 模块读取 `pressure_neumann` 的接口。它只跳过数学上与 `regular` 完全相同的重复构造。

优化后重新运行同一静水基准算例，原始耗时输出为：

第 1 步：

```text
[timing] step=1 time=0.001 neighbor=0.161391 time_control=0.0156142 free_surface=2.50125 lsmps=5.58869 provisional=1.28565 ppe=6.83222 correction=0.445202 diagnostics=0.0157105 vtk=0.0686648 total=16.9144
```

第 2 步：

```text
[timing] step=2 time=0.002 neighbor=0.163901 time_control=0.0154702 free_surface=2.53768 lsmps=5.70405 provisional=1.30169 ppe=6.88736 correction=0.452807 diagnostics=0.0158247 vtk=0.0758985 total=17.1547
```

本轮优化后平均耗时统计：

| 模块 | 平均耗时 s | 占总耗时 |
|---|---:|---:|
| PPE 构造与求解 | 6.85979 | 40.27% |
| LSMPS 矩阵构造 | 5.64637 | 33.15% |
| 自由面识别 | 2.51947 | 14.79% |
| Provisional | 1.29367 | 7.59% |
| Correction | 0.44900 | 2.64% |
| 邻域搜索 | 0.16265 | 0.95% |
| VTK 输出 | 0.07228 | 0.42% |
| 时间步控制 | 0.01554 | 0.09% |
| 诊断统计 | 0.01577 | 0.09% |
| 总耗时 | 17.03455 | 100% |

与上一轮 PPE 散度优化后的基准相比：

- LSMPS 平均耗时：`7.17414 s` -> `5.64637 s`，下降约 `21.30%`。
- 总单步耗时：`18.64495 s` -> `17.03455 s`，下降约 `8.64%`。

与最初基准相比：

- 总单步耗时：`21.17425 s` -> `17.03455 s`，累计下降约 `19.55%`。

### 5.2 已完成优化：LSMPS 诊断配置开关

优化时间：2026-06-23

优化位置：

- `src/core/simulation_config.hpp`
- `src/config/config_reader.cpp`
- `src/lsmps/lsmps_matrices.cpp`
- `cases/hydrostatic_box/config.ini`

新增配置项：

```ini
[lsmps]
diagnostics_enabled = false
```

当 `diagnostics_enabled = true` 时，保持原有完整诊断路径：对 raw moment matrix 做特征值分解、秩判断、条件数估计，并根据阈值给出 `Valid`、`IllConditioned` 或失败状态。

当 `diagnostics_enabled = false` 时，跳过特征值分解和条件数估计，只使用 LU 可逆性判断矩阵是否可用，然后计算逆矩阵。该模式不改变矩阵构造公式和后续算子离散，但会减少坏粒子的详细诊断信息，VTK 中条件数等字段不再具有诊断意义。

静水基准算例使用该快速模式。

## 6. 自由面识别瓶颈分析

相关代码：

- `src/free_surface/free_surface_detector.cpp`

当前自由面识别平均耗时约 2.51 s，占总时间约 11.8%。主要原因包括：

- 对所有流体粒子执行 cubed sphere 方向采样。
- 每个粒子会遍历方向、流体邻居和壁面邻居来构造遮挡/开口信息。
- 静水算例中自由面变化很小，但当前每步全量重算。

优先优化建议：

1. 已完成：增加保守数密度粗筛，先排除明显内部粒子，再执行完整 `R_open + R_cone`。
2. 降低非关键算例中的 `cubed_sphere_q`，或让其按精度等级配置。
3. 对静水或小变形算例，允许低频更新自由面状态。
4. 对粒子循环增加 OpenMP 并行。

### 6.1 已完成优化：自由面识别数密度粗筛

优化时间：2026-06-23

优化位置：

- `src/free_surface/free_surface_detector.cpp`
- `cases/hydrostatic_box/config.ini`

粗筛逻辑：

1. 对无壁面邻居的流体粒子计算加权流体数密度。
2. 用当前粒子系统中的最大无壁面数密度作为内部参考数密度。
3. 当粒子同时满足以下条件时，直接判定为内部粒子并跳过 `R_open + R_cone`：
   - 无壁面邻居。
   - 流体邻居数大于 `splash_max_fluid_neighbors`。
   - 加权数密度不小于 `number_density_ratio * reference_number_density`。

本次静水算例设置：

```ini
[free_surface]
number_density_ratio = 0.98
```

该阈值比此前默认的 `0.90` 更保守，目的是避免把真正自由面粒子误筛为内部粒子。测试中曾发现仅用数密度会误筛 splash 粒子，因此最终粗筛额外要求流体邻居数大于 `splash_max_fluid_neighbors`。

### 6.2 已完成优化：PPE 按行批量装配

优化时间：2026-06-23

优化位置：

- `src/pressure_poisson/pressure_poisson.cpp`

优化前，PPE 矩阵装配对每个非零系数单独调用一次 `MatSetValues`。优化后，每个粒子行先收集 `columns` 和 `values`，再对该行执行一次批量 `MatSetValues`。

该优化不改变 PPE 矩阵数值和右端项，只减少 PETSc 装配阶段的小粒度函数调用。当前 8000 流体粒子静水基准中 PPE 总耗时未明显下降，说明该算例下 PPE 时间主要仍由求解器迭代和矩阵整体装配/预条件过程控制。

本轮三项优化后重新运行同一静水基准算例，原始耗时输出为：

第 1 步：

```text
[timing] step=1 time=0.001 neighbor=0.163212 time_control=0.0155237 free_surface=1.70161 lsmps=4.51146 provisional=1.29959 ppe=6.94685 correction=0.452294 diagnostics=0.0157623 vtk=0.080858 total=15.1872
```

第 2 步：

```text
[timing] step=2 time=0.002 neighbor=0.161024 time_control=0.0152836 free_surface=1.72168 lsmps=4.58467 provisional=1.31661 ppe=6.99339 correction=0.454852 diagnostics=0.0154167 vtk=0.0741251 total=15.337
```

本轮优化后平均耗时统计：

| 模块 | 平均耗时 s | 占总耗时 |
|---|---:|---:|
| PPE 构造与求解 | 6.97012 | 45.67% |
| LSMPS 矩阵构造 | 4.54806 | 29.80% |
| 自由面识别 | 1.71165 | 11.22% |
| Provisional | 1.30810 | 8.57% |
| Correction | 0.45357 | 2.97% |
| 邻域搜索 | 0.16212 | 1.06% |
| VTK 输出 | 0.07749 | 0.51% |
| 时间步控制 | 0.01540 | 0.10% |
| 诊断统计 | 0.01559 | 0.10% |
| 总耗时 | 15.26210 | 100% |

与上一轮优化后相比：

- LSMPS 平均耗时：`5.64637 s` -> `4.54806 s`，下降约 `19.45%`。
- 自由面识别平均耗时：`2.51947 s` -> `1.71165 s`，下降约 `32.06%`。
- 总单步耗时：`17.03455 s` -> `15.26210 s`，下降约 `10.41%`。

与最初基准相比：

- 总单步耗时：`21.17425 s` -> `15.26210 s`，累计下降约 `27.92%`。

## 7. 建议优化路线

第一阶段，低风险优化：

1. 已完成：优化 PPE `computeDivergence` 的重复数组构造。
2. 已完成：对无壁面邻居粒子复用 `regular` 作为 `pressure_neumann`。
3. 不采用跨时间步 LSMPS 矩阵复用；为保持离散精度，每个时间步仍重新计算 LSMPS 矩阵。

第二阶段，并行和缓存：

1. 暂不实现 OpenMP 粒子并行。
2. 暂不实现邻居几何 stencil 缓存，以控制内存占用。
3. 后续若继续优化 CPU 性能，应优先评估 PETSc 求解器和预条件器参数。

第三阶段，PPE 求解器优化：

1. 测试 PETSc 不同 KSP/PC 组合。
2. 引入多进程 PETSc 装配和求解。
3. 在矩阵结构变化较小时复用 PETSc Mat 预分配结构。

## 8. 当前结论

本次约 8000 流体粒子的静水基准表明，经过三轮优化后，当前平均单步耗时约 `15.26210 s`。相比最初基准 `21.17425 s`，总单步耗时累计下降约 `27.92%`。

最新热点排序为：

1. PPE 构造与求解：`6.97012 s`，约 `45.67%`。
2. LSMPS 矩阵构造：`4.54806 s`，约 `29.80%`。
3. 自由面识别：`1.71165 s`，约 `11.22%`。
4. Provisional：`1.30810 s`，约 `8.57%`。

LSMPS 和自由面识别已经获得较明显加速。PPE 仍是当前最大瓶颈，且 PPE 按行批量装配对该算例收益有限，说明后续 PPE 优化重点应转向 PETSc KSP/PC 参数、预条件器和求解器路径。

后续加速不应直接从 VTK 或邻域搜索开始。更合理的顺序是：

1. PETSc KSP/PC 参数调优。
2. 评估 PPE 矩阵装配与预分配策略。
3. 继续在不复用跨步 LSMPS 矩阵的前提下优化单步 LSMPS 构造。
4. 保守调整自由面粗筛阈值，并持续通过 VTK 验证避免漏判自由面。
