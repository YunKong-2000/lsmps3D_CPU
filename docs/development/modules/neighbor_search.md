# 邻居粒子搜索模块实现说明

本文档记录当前三维邻域搜索模块的实现细节、接口约定和测试覆盖。模块源码位于 `src/neighbor/`。

## 1. 目标与边界

邻域搜索模块用于在给定支持域半径 `support_radius` 下，为每个粒子找出距离不超过该半径的其他粒子。模块只依赖 `ParticleSet` 和几何位置，不依赖具体物理方程。

当前实现支持：

- 流体粒子和壁面粒子共同参与搜索。
- 为每个粒子生成邻居索引列表。
- 分别保存总邻居索引、流体邻居索引和壁面邻居索引。
- 将总邻居数量回写到 `ParticleSet::neighborCounts()`，供 VTK 输出、自由面识别和诊断使用。
- 将流体邻居数量和壁面邻居数量分别回写到 `fluidNeighborCounts()` 和 `wallNeighborCounts()`。
- 负坐标粒子的 cell 映射。

当前实现暂不包含：

- 周期边界。
- 动态裁剪的规则三维数组网格。
- 并行构建。
- 按粒子类型过滤邻居。

## 2. 模块结构

### 2.1 `CellLinkedList`

文件：

- `src/neighbor/cell_linked_list.hpp`
- `src/neighbor/cell_linked_list.cpp`

职责：

- 保存背景网格原点 `origin` 和 cell 尺寸 `cell_size`。
- 将粒子坐标映射到整数三维 cell 索引。
- 使用哈希表保存非空 cell 到粒子索引列表的映射。
- 提供 `particlesInCell()` 查询单个 cell 内的候选粒子。

当前 cell 尺寸取 `support_radius`。这样任意有效邻居只可能位于目标粒子所在 cell 及周围 26 个相邻 cell 中。

### 2.2 `NeighborSearch`

文件：

- `src/neighbor/neighbor_search.hpp`
- `src/neighbor/neighbor_search.cpp`

职责：

- 根据 `support_radius` 构建 `CellLinkedList`。
- 对每个粒子访问 27 个候选 cell。
- 对候选粒子做真实欧氏距离筛选。
- 返回 `TypedNeighborList`，其中包含：
  - `fluid`：流体邻居索引。
  - `wall`：壁面邻居索引。
- 兼容保留 `NeighborList`，类型为 `std::vector<std::vector<std::size_t>>`，由分类邻居列表按需合并得到。
- 可通过 `countNeighborsByType()` 统计总邻居、流体邻居和壁面邻居数量。
- 可将分类邻居数量写回 `ParticleSet`。

## 3. Cell 映射

粒子位置 `position` 到 cell 索引的映射为：

```text
cell_x = floor((position.x - origin.x) / cell_size)
cell_y = floor((position.y - origin.y) / cell_size)
cell_z = floor((position.z - origin.z) / cell_size)
```

使用 `floor` 而不是截断，可以正确处理负坐标。例如当 `origin.x = 0` 且 `cell_size = 1` 时：

- `x = 0.0` 映射到 cell `0`
- `x = 1.0` 映射到 cell `1`
- `x = -0.1` 映射到 cell `-1`

哈希表只保存非空 cell，因此不需要预先知道全局网格尺寸，也不需要给大空域分配内存。

## 4. 搜索流程

对每次邻域搜索：

1. 清空 cell 哈希表。
2. 遍历所有粒子，将粒子索引插入对应 cell。
3. 对每个目标粒子 `i`：
   - 计算目标粒子所在 cell。
   - 遍历 `dx, dy, dz in [-1, 0, 1]` 的 27 个候选 cell。
   - 跳过自身粒子。
   - 计算候选粒子和目标粒子的距离平方。
   - 若 `distance_squared <= support_radius^2`，加入总邻居列表。
   - 根据候选粒子类型，同时加入流体邻居列表或壁面邻居列表。
4. 对每个粒子的邻居索引排序，保证结果稳定，方便测试和调试。

距离判定包含半径边界，即距离正好等于 `support_radius` 的粒子会被视为邻居。

## 5. 复杂度

设粒子数为 `N`，每个 cell 内平均候选粒子数较小且近似稳定，则：

- 构建 cell linked-list：`O(N)`
- 查询邻居：约 `O(N * k)`，其中 `k` 是 27 个候选 cell 内的平均候选数量。

在粒子分布近似均匀、cell 尺寸取支持半径时，`k` 通常远小于 `N`。最坏情况下，如果所有粒子都落入同一区域，复杂度仍会退化到 `O(N^2)`。

## 6. 与其他模块的关系

- `ParticleSet`：邻域搜索读取粒子位置，并将邻居数量写入 `neighborCounts()`、`fluidNeighborCounts()` 和 `wallNeighborCounts()`。
- `VtkWriter`：默认输出 `neighbor_count`、`fluid_neighbor_count` 和 `wall_neighbor_count` 字段，可用于 ParaView 检查搜索结果。
- 自由面识别：后续可直接使用 `TypedNeighborList::fluid`、`TypedNeighborList::wall` 或分类邻居数量进行阈值判定。
- LSMPS 算子：后续将使用 `TypedNeighborList` 遍历局部支持域，构造 type-A 矩矩阵和差分算子。

## 7. 测试覆盖

测试文件为 `tests/neighbor_search_test.cpp`，当前覆盖：

- 3x3x3 规则粒子分布。
- 与暴力 `O(N^2)` 搜索结果逐粒子对比。
- 角点邻居数量为 3。
- 中心点邻居数量为 6。
- 邻居数量回写到 `ParticleSet::neighborCounts()`。
- 流体粒子和壁面粒子共同参与搜索。
- 小型混合测试直接检查流体邻居索引列表和壁面邻居索引列表。
- 343 个粒子的较大混合测试：5x5x5 流体块加一层六面壁面粒子壳。
- 较大混合测试中，中心流体粒子的流体邻居数为 6、壁面邻居数为 0。
- 较大混合测试中，靠壁流体粒子的流体邻居数为 5、壁面邻居数为 1。
- 较大混合测试直接检查分类邻居索引列表大小。
- 总邻居数量等于流体邻居数量与壁面邻居数量之和。
- 50x50x50 规则流体粒子大规模测试，共 125000 个流体粒子。
- 大规模测试在流体块外添加单层壁面粒子壳，共 15608 个壁面粒子。
- 大规模测试支持域半径取 `3.1 * particle_spacing`。
- 大规模测试将所有粒子的 `neighbor_count`、`fluid_neighbor_count` 和 `wall_neighbor_count` 输出到 `output/neighbor_search_50x50x50.vtk`，用于 ParaView 可视化验收。
- 大规模测试通过局部整数偏移枚举逐流体粒子校验解析流体邻居数和壁面邻居数。
- 远距离粒子不会被错误纳入邻居。
- 负坐标、边界坐标的 cell 映射。
- 非法支持半径会抛出异常。

## 8. 后续扩展点

- 根据后续算法需要，为邻居列表缓存相对位移和距离平方。
- 增加邻居距离、相对位移缓存，减少 LSMPS 阶段重复计算。
- 支持周期边界或开放边界 ghost particle 策略。
- 对大规模算例加入并行构建和查询。
- 将邻居数量、最小/最大/平均邻居数写入 `SimulationState` 诊断字段。
