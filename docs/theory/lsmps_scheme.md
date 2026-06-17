# lsmps格式推导
## 泰勒展开式
对于物理量$\phi$和观测流体粒子i,以及邻域中流体粒子j,他们的物理量通过泰勒展开式表达为,
$$
\phi_j = \phi_i + \sum_{k=1}^{n} \frac{1}{k!}(\mathbf{r}_{ij} \cdot \nabla)^k \phi_i
$$
其中$\mathbf{r}_{ij}$为i,j粒子的距离矢量。上述泰勒展开式的展开形式为
$$
\phi_{j} = \phi_i + \frac{\partial \phi}{\partial x}x_{ij} + \frac{\partial \phi}{\partial y}y_{ij} + \frac{\partial \phi}{\partial z}z_{ij} + \frac{1}{2}\frac{\partial^2 \phi}{\partial x^2}x_{ij}^2 + \frac{1}{2}\frac{\partial^2 \phi}{\partial y^2}y_{ij}^2 + \frac{1}{2}\frac{\partial^2 \phi}{\partial z^2}z_{ij}^2 + \frac{\partial^2 \phi}{\partial x \partial y}x_{ij}y_{ij} + \frac{\partial^2 \phi}{\partial x \partial z}x_{ij}z _{ij} + \frac{\partial^2 \phi}{\partial y \partial z}y_{ij}z_{ij} + O(|r_{ij}|^3)
$$
## type-A格式
设计一个归一化的基向量$\mathbf{p}_{ij}$,
$$
\mathbf{p}_{ij}=[\frac{x_{ij}}{r_e}, \frac{y_{ij}}{r_e}, \frac{z_{ij}}{r_e}, \frac{x_{ij}^2}{2r_e^2},\frac{y_{ij}^2}{2r_e^2}, \frac{z_{ij}^2}{2r_e^2}, \frac{x_{ij}y_{ij}}{r_e^2}, \frac{x_{ij}z_{ij}}{r_e^2}, \frac{y_{ij}z_{ij}}{r_e^2}]^T
$$
对应设计一个导数特征向量$\mathbf{D}_i$,
$$
\mathbf{D}_i=[r_e\frac{\partial \phi}{\partial x}, r_e\frac{\partial \phi}{\partial y}, r_e\frac{\partial \phi}{\partial z}, r_e^2\frac{\partial^2 \phi}{\partial x^2}, r_e^2\frac{\partial^2 \phi}{\partial y^2}, r_e^2\frac{\partial^2 \phi}{\partial z^2},r_e^2\frac{\partial^2 \phi}{\partial x \partial y}, r_e^2\frac{\partial^2 \phi}{\partial x \partial z}, r_e^2\frac{\partial^2 \phi}{\partial y \partial z}]_i^T
$$
于是可以获得一个近似的对物理量差值$\phi_j - \phi_i$的近似估计，
$$
\phi_j - \phi_i \approx \mathbf{p}_{ij}^T\mathbf{D}_i
$$

于是可以对次近似估计建立误差目标函数$J(\mathbf{D}_i)$
$$
J(\mathbf{D}_i)=\sum_{j \neq i}w_{ij}[\mathbf{p}_{ij}^T\mathbf{D}_i - (\phi_j - \phi_i)]^2
$$
求极值问题$\frac{\partial J(\mathbf{D}_i)}{\partial \mathbf{D}_i}=0$,则可以得到
$$
\frac{\partial J(\mathbf{D}_i)}{\partial \mathbf{D}_i}=\sum_{j \neq i}w_{ij}2\mathbf{p}_{ij}[\mathbf{p}_{ij}^T\mathbf{D}_i - (\phi_j - \phi_i)]=0
$$
得到
$$
\sum_{j \neq i}w_{ij}\mathbf{p}_{ij}\mathbf{p}_{ij}^T\mathbf{D}_i=\sum_{j \neq i}w_{ij}\mathbf{p}_{ij}(\phi_j - \phi_i) \\
\mathbf{D}_i=(\sum_{j \neq i}w_{ij}\mathbf{p}_{ij}\mathbf{p}_{ij}^T)^{-1}\sum_{j \neq i}w_{ij}\mathbf{p}_{ij}(\phi_j - \phi_i)
$$
于是可以得到局部下使得近似误差最小的一组导数，
令$M_i=(\sum_{j \neq i}w_{ij}\mathbf{p}_{ij}\mathbf{p}_{ij}^T)^{-1}$,也即是i粒子对应的矩矩阵。可以看出矩矩阵的计算仅与i粒子周围的粒子分布位置相关，与具体物理量无关。
则可以得到lsmps typeA格式下的各个导数的离散形式，
$$
\frac{\partial \phi_i}{\partial x}=\frac{1}{r_e}\mathbf{D}_{i,0}=\frac{1}{r_e}\sum_{j \neq i}w_{ij}M_{i,0}\mathbf{p}_{ij}(\phi_j - \phi_i) \\
\frac{\partial \phi_i}{\partial y}=\frac{1}{r_e}\mathbf{D}_{i,1}=\frac{1}{r_e}\sum_{j \neq i}w_{ij}M_{i,1}\mathbf{p}_{ij}(\phi_j - \phi_i) \\
\frac{\partial \phi_i}{\partial z}=\frac{1}{r_e}\mathbf{D}_{i,2}=\frac{1}{r_e}\sum_{j \neq i}w_{ij}M_{i,2}\mathbf{p}_{ij}(\phi_j - \phi_i) \\
\frac{\partial^2 \phi_i}{\partial x^2}=\frac{1}{r_e^2}\mathbf{D}_{i,3}=\frac{1}{r_e^2}\sum_{j \neq i}w_{ij}M_{i,3}\mathbf{p}_{ij}(\phi_j - \phi_i) \\
\frac{\partial^2 \phi_i}{\partial y^2}=\frac{1}{r_e^2}\mathbf{D}_{i,4}=\frac{1}{r_e^2}\sum_{j \neq i}w_{ij}M_{i,4}\mathbf{p}_{ij}(\phi_j - \phi_i) \\
\frac{\partial^2 \phi_i}{\partial z^2}=\frac{1}{r_e^2}\mathbf{D}_{i,5}=\frac{1}{r_e^2}\sum_{j \neq i}w_{ij}M_{i,5}\mathbf{p}_{ij}(\phi_j - \phi_i) \\
$$
可以看出不同的偏导数计算对应了矩矩阵中不同的行。
自然可以得到不同算子的离散形式，
$$
<\nabla \phi>_i=\frac{1}{r_e}\sum_{j \neq i}w_{ij}
\begin{bmatrix}
M_{i,0} \\
M_{i,1} \\
M_{i,2} \\
\end{bmatrix} \mathbf{p}_{ij}(\phi_j - \phi_i) \\[5pt]
<\nabla \cdot \mathbf{\phi}>_i=\frac{1}{r_e}\sum_{j \neq i}w_{ij}
\begin{bmatrix}
M_{i,0} \\
M_{i,1} \\
M_{i,2} \\
\end{bmatrix} \mathbf{p}_{ij} \cdot (\mathbf{\phi}_j - \mathbf{\phi}_i) \\[5pt]
<\nabla^2 \phi>_i=\frac{1}{r_e^2}\sum_{j \neq i}w_{ij}(M_{i,3}+M_{i,4}+M_{i,5})\mathbf{p}_{ij}(\phi_j - \phi_i)
$$

## type-B格式
type-B格式用于说明另一种lsmps拟合思路。当前程序在PPE构造中不采用type-B格式，本节只作为数学背景保留。

对于泰勒展开式
$$
\phi_j = \phi_i + \sum_{k=1}^{n} \frac{1}{k!}(\mathbf{r}_{ij} \cdot \nabla)^k
$$
也可以写成，
$$
\phi_j - \phi_i = \sum_{k=1}^{n} \frac{1}{k!}(\mathbf{r}_{ij} \cdot \nabla)^k
$$
区别在于前者将i粒子的物理量$\phi_i$视为一个已知确定的值，type-A格式本质是在这个前提下在局部进行拟合。但对于i粒子物理量本身就不确定待求解的情况，就可以将i粒子物理量本身作为拟合的一个目标。也就引出了后者type-B格式。
这时候修改基向量$\mathbf{p}_{ij}$,
$$
\mathbf{p}_{ij}=[1,\frac{x_{ij}}{r_e}, \frac{y_{ij}}{r_e}, \frac{z_{ij}}{r_e}, \frac{x_{ij}^2}{2r_e^2},\frac{y_{ij}^2}{2r_e^2}, \frac{z_{ij}^2}{2r_e^2}, \frac{x_{ij}y_{ij}}{r_e^2}, \frac{x_{ij}z_{ij}}{r_e^2}, \frac{y_{ij}z_{ij}}{r_e^2}]^T
$$
对应的导数特征向量$\mathbf{D}_i$为，
$$
\mathbf{D}_i=[\phi, r_e\frac{\partial \phi}{\partial x}, r_e\frac{\partial \phi}{\partial y}, r_e\frac{\partial \phi}{\partial z}, r_e^2\frac{\partial^2 \phi}{\partial x^2}, r_e^2\frac{\partial^2 \phi}{\partial y^2}, r_e^2\frac{\partial^2 \phi}{\partial z^2},r_e^2\frac{\partial^2 \phi}{\partial x \partial y}, r_e^2\frac{\partial^2 \phi}{\partial x \partial z}, r_e^2\frac{\partial^2 \phi}{\partial y \partial z}]_i^T
$$
此时目标函数$J(\mathbf{D}_i)$也有相应变化，
$$
J(\mathbf{D}_i)=\sum_{j}w_{ij}[\mathbf{p}_{ij}^T\mathbf{D}_i - \phi_j]^2
$$
可以看到，此时由于$\phi_i$不再是已知量，所以目标函数中仅有$\phi_j$。此时$M_i=(\sum_{j} w_{ij}\mathbf{p}_{ij}\mathbf{p}_{ij}^T)^{-1}$而极值对应的导数特征向量也变为，
$$
\mathbf{D}_{i}=\sum_{j}w_{ij}M_i\mathbf{p}_{ij}\phi_j
$$
可以看到，求和符号的下缀不再约束$j \neq i$,也就是说在type-B格式下，i粒子本身将会以加权的形式参与到拟合中，显然i粒子对于增加物理量的权重是最大的。后续对算子进行离散时也出现相应的变化。
$$
<\nabla \phi>_i=\frac{1}{r_e}\sum_{j}w_{ij}
\begin{bmatrix}
M_{i,1} \\
M_{i,2} \\
M_{i,3} \\
\end{bmatrix} \mathbf{p}_{ij}\phi_j \\[5pt]
<\nabla \cdot \mathbf{\phi}>_i=\frac{1}{r_e}\sum_{j}w_{ij}
\begin{bmatrix}
M_{i,1} \\
M_{i,2} \\
M_{i,3} \\
\end{bmatrix} \mathbf{p}_{ij} \cdot \mathbf{\phi}_j\\[5pt]
<\nabla^2 \phi>_i=\frac{1}{r_e^2}\sum_{j}w_{ij}(M_{i,4}+M_{i,5}+M_{i,6})\mathbf{p}_{ij}\phi_j
$$

## 边界条件额外约束
lsmps方法中，部分边界条件可以通过在目标函数中添加额外的约束来实现。这种实现是局部的，也是更加直接的。

当前程序中需要区分两类处理方式：

- 第一类边界条件，例如自由面压力$P=P_{\mathrm{fs}}$，不通过修改lsmps算子或矩矩阵实现，而是在最终线性方程组中直接替换对应粒子的方程行。
- 第二类边界条件，例如壁面压力法向导数$\nabla P \cdot \mathbf{n}=C$，可以作为额外约束加入局部最小二乘目标函数，从而影响对应粒子的矩矩阵和算子系数。

因此，后续PPE构造中不会为了自由面压力第一类边界条件引入额外的$\lambda_{FreeSurface}\mathbf{c}\mathbf{c}^T$项，也不会为自由面粒子使用type-B格式构造压力拉普拉斯算子。

### 第一类边界条件

第一类边界条件直接指定物理量本身：

$$
\phi_i = \phi_D
$$

对于PPE中的自由面压力边界条件，通常取

$$
P_i = P_{\mathrm{fs}} = 0
$$

程序实现时，应在压力线性系统装配阶段将自由面粒子的方程替换为

$$
A_{ii}=1,\qquad A_{ij}=0\ (j\neq i),\qquad b_i=P_{\mathrm{fs}}
$$

这种处理不改变lsmps局部算子本身。对于其他粒子的压力泊松方程，如果邻域中包含自由面粒子，可以保留自由面压力未知量并依靠其自身的Dirichlet方程约束，也可以在装配时将已知的$P_{\mathrm{fs}}$移到右端项。初始程序建议采用前者，即所有流体粒子都保留压力自由度，自由面粒子使用Dirichlet方程行。

### 第二类边界条件
施加第二类边界条件$\frac{\partial \phi}{\partial \mathbf{n}}=\nabla \phi \cdot \mathbf{n}=C$的方法应该类似：使用one-hot向量或者变形后的one-hot向量将$\mathbf{D}_i$中对应的偏导数提取出来然后作为额外的约束加入目标函数。但是实际上的$\mathbf{n}$方向通常不与$x,y,z$方向重合，所以无法直接提取出对应的方向导数，需要先对泰勒展开式做一定变化，将其转到$\mathbf{n}$方向后再进行离散。  
在局部空间中定义的导数特征向量$\mathbf{D}_i$的值是恒定的常数，对于type-A格式，i粒子的物理量也是恒定的
所以
$$
\nabla \phi_j=\nabla (\phi_i + \mathbf{p}_{ij}^T\mathbf{D}_i)=(\nabla \mathbf{p}_{ij})^T\mathbf{D}_i
$$
type-A格式下的基向量$\mathbf{p}_{ij}=[\frac{x_{ij}}{r_e}, \frac{y_{ij}}{r_e}, \frac{z_{ij}}{r_e}, \frac{x_{ij}^2}{2r_e^2},\frac{y_{ij}^2}{2r_e^2}, \frac{z_{ij}^2}{2r_e^2}, \frac{x_{ij}y_{ij}}{r_e^2}, \frac{x_{ij}z_{ij}}{r_e^2}, \frac{y_{ij}z_{ij}}{r_e^2}]^T$
则，
$$
\nabla \mathbf{p}_{ij} = 
\begin{bmatrix}
\frac{1}{r_e},0,0,\frac{x_{ij}}{r_e^2},0,0,\frac{y_{ij}}{r_e^2},\frac{z_{ij}}{r_e^2},0 \\[5pt]
0,\frac{1}{r_e},0,0,\frac{y_{ij}}{r_e^2},0,\frac{x_{ij}}{r_e^2},0,\frac{z_{ij}}{r_e^2} \\[5pt]
0,0,\frac{1}{r_e},0,0,\frac{z_{ij}}{r_e^2},0,\frac{x_{ij}}{r_e^2},\frac{y_{ij}}{r_e^2}
\end{bmatrix}^T
$$
进一步，
$$
\nabla \mathbf{p}_{ij} \cdot \mathbf{n}= \nabla \mathbf{p}_{ij} \cdot [n_x,n_y,n_z]=
[\frac{n_x}{r_e},\frac{n_y}{r_e},\frac{n_z}{r_e},
\frac{n_xx_{ij}}{r_e^2},\frac{n_yy_{ij}}{r_e^2},\frac{n_zz_{ij}}{r_e^2},
\frac{n_xy_{ij}+n_yx_{ij}}{r_e^2},
\frac{n_xz_{ij}+n_zx_{ij}}{r_e^2},
\frac{n_yz_{ij}+n_zy_{ij}}{r_e^2}
]^T
$$
令新的基向量
$$\mathbf{q}_{ij}=r_e\nabla \mathbf{p}_{ij} \cdot \mathbf{n}
=[n_x, n_y, n_z, \frac{n_xx_{ij}}{r_e},\frac{n_yy_{ij}}{r_e},\frac{n_zz_{ij}}{r_e},
\frac{n_xy_{ij}+n_yx_{ij}}{r_e},
\frac{n_xz_{ij}+n_zx_{ij}}{r_e},
\frac{n_yz_{ij}+n_zy_{ij}}{r_e}
]^T
$$
将原本的边界条件变形为$r_e\nabla \phi \cdot \mathbf{n}=r_eC$，
则$r_e \nabla \phi \cdot \mathbf{n}=\mathbf{q}_{ij}^T\mathbf{D}_i=r_eC$。
现在将边界条件约束其加入至目标函数，
$$
J(\mathbf{D}_i)=\sum_{j \in Fluid,\ j \neq i}w_{ij}[\mathbf{p}_{ij}^T\mathbf{D}_i-(\phi_j - \phi_i)]^2+\sum_{j \in Wall}w_{ij}(\mathbf{q}_{ij}^T\mathbf{D}_i-r_eC)^2
$$
求极值问题，
$$
\frac{\partial J(\mathbf{D}_i)}{\partial \mathbf{D}_i}=\sum_{j \in Fluid,\ j \neq i}w_{ij}\mathbf{p}_{ij}(\mathbf{p}_{ij}^T\mathbf{D}_{i}-(\phi_j - \phi_i))+\sum_{j \in Wall}w_{ij}\mathbf{q}_{ij}(\mathbf{q}_{ij}^T\mathbf{D}_i-r_eC)=0
$$
最终得到
$$
\mathbf{D}_i=(\sum_{j \in Fluid,\ j \neq i}w_{ij}\mathbf{p}_{ij}\mathbf{p}_{ij}^T+\sum_{j \in Wall}w_{ij}\mathbf{q}_{ij}\mathbf{q}_{ij}^T)^{-1}(\sum_{j \in Fluid,\ j \neq i}w_{ij}\mathbf{p}_{ij}(\phi_j - \phi_i)+\sum_{j \in Wall}w_{ij}\mathbf{q}_{ij}r_eC)
$$
令$M=(\sum_{j \in Fluid,\ j \neq i}w_{ij}\mathbf{p}_{ij}\mathbf{p}_{ij}^T+\sum_{j \in Wall}w_{ij}\mathbf{q}_{ij}\mathbf{q}_{ij}^T)^{-1}$，则$M$为代入了第二类边界条件的矩矩阵。
以近壁面流体粒子的压力第二类边界为例，离散的算子形式为，
$$
<\nabla P>_i=\frac{1}{r_e}(\sum_{j \in Fluid,\ j \neq i}w_{ij}
\begin{bmatrix}
M_{i,0} \\
M_{i,1} \\
M_{i,2} \\
\end{bmatrix} \mathbf{p}_{ij}(P_j - P_i) 
+ \sum_{j \in Wall}w_{ij}
\begin{bmatrix}
M_{i,0} \\
M_{i,1} \\
M_{i,2} \\
\end{bmatrix} \mathbf{q}_{ij}r_eC)
\\[5pt]
<\nabla^2 P>_i=\frac{1}{r_e^2}(\sum_{j \in Fluid,\ j \neq i}w_{ij}(M_{i,3}+M_{i,4}+M_{i,5})\mathbf{p}_{ij}(P_j - P_i) + \sum_{j \in Wall}w_{ij}(M_{i,3}+M_{i,4}+M_{i,5})\mathbf{q}_{ij}r_eC)
$$

type-B格式也可以推导类似的第二类边界条件约束形式，但当前程序的PPE构造不采用type-B格式。后续实现PPE时，应只使用上面的type-A边界约束形式；type-B部分仅作为lsmps数学格式的背景说明，不作为压力方程装配依据。


## 压力泊松方程(PPE)离散
lsmps中PPE的系数构造和右边源项的构建同时涉及速度散度和压力拉普拉斯算子的离散，且不同边界条件会影响线性系统的装配方式。以下部分详细说明PPE的系数矩阵和右边源项的具体计算方法。
### 压力泊松方程原始形式
PPE的原始形式为，
$$
\frac{1}{\rho}\nabla^2 P_{i}^{k+1}=\frac{1}{\Delta t}\nabla \cdot \mathbf{u}_i^*
$$
其中$P_i^{k+1}$为待求解的i粒子的压力大小，$\mathbf{u}_i^*$为i粒子的临时速度。
### 临时速度散度的离散
临时速度的散度是右边源项的主要构成部分。在自由面流动模拟中，壁面处通常使用不可滑移边界，也即壁面粒子和靠近壁面粒子的流体粒子的相对速度为零。速度边界条件不通过修改lsmps散度算子实现，实际计算时直接将壁面粒子的临时速度代入type-A格式的常规散度离散即可。不过需要注意的是，壁面粒子的速度也必须是临时速度，而不是真实速度，壁面粒子的临时速度计算方法为，
$$
\mathbf{u}_{Wall}^*= \mathbf{u}_{Wall}+\mathbf{g}\Delta t
$$
也即是，壁面粒子在计算临时速度需要考虑重力的作用。
最终临时速度的离散形式为，
$$
<\nabla \cdot \mathbf{u}^*>_i=\frac{1}{r_e}\sum_{j \neq i}w_{ij}
\begin{bmatrix}
M_{i,0} \\
M_{i,1} \\
M_{i,2} \\
\end{bmatrix} \mathbf{p}_{ij} \cdot (\mathbf{u}_j^*-\mathbf{u}_i^*)
$$
其中$M_i=(\sum_{j \neq i} w_{ij}\mathbf{p}_{ij} \mathbf{p}_{ij}^T)^{-1}$。可见在速度散度的离散中，壁面粒子基本视为了一种特殊的流体粒子，除了其临时速度的计算方法和流体粒子略有区别外和流体粒子基本相同。

### 压力拉普拉斯算子的离散
前文已经讨论过边界条件的实现。当前程序中，PPE的压力拉普拉斯算子统一使用type-A格式离散，不再使用type-B格式。

对非自由面的流体粒子，压力拉普拉斯算子的基本形式为

$$
<\nabla^2 P>_i=
\frac{1}{r_e^2}
\sum_{j \in Fluid,\ j \neq i}
w_{ij}
(M_{i,3}+M_{i,4}+M_{i,5})
\mathbf{p}_{ij}
(P_j-P_i)
$$

其中，矩矩阵$M_i$根据该粒子是否需要施加壁面压力第二类边界条件分为两种情况。

不包含壁面第二类边界条件时：

$$
M_i =
\left(
\sum_{j \in Fluid,\ j \neq i}
w_{ij}\mathbf{p}_{ij}\mathbf{p}_{ij}^T
\right)^{-1}
$$

包含壁面第二类边界条件时：

$$
M_i =
\left(
\sum_{j \in Fluid,\ j \neq i}
w_{ij}\mathbf{p}_{ij}\mathbf{p}_{ij}^T
+
\sum_{j \in Wall}
w_{ij}\mathbf{q}_{ij}\mathbf{q}_{ij}^T
\right)^{-1}
$$

此时压力拉普拉斯算子写为

$$
<\nabla^2 P>_i=
\frac{1}{r_e^2}
\sum_{j \in Fluid,\ j \neq i}
w_{ij}
(M_{i,3}+M_{i,4}+M_{i,5})
\mathbf{p}_{ij}
(P_j-P_i)
+
\frac{1}{r_e^2}
\sum_{j \in Wall}
w_{ij}
(M_{i,3}+M_{i,4}+M_{i,5})
\mathbf{q}_{ij}
r_e C_j
$$

对压力边界，$C_j$通常取

$$
C_j = \rho \mathbf{g}\cdot \mathbf{n}_j
$$

其中$\mathbf{n}_j$为壁面粒子法向，方向应与程序中的边界约定保持一致。

自由面压力第一类边界条件$P=P_{\mathrm{fs}}$不进入上述矩矩阵构造，也不改变lsmps压力拉普拉斯算子的离散格式。自由面粒子的PPE方程在装配线性系统时直接替换为Dirichlet方程行。

### PPE系数矩阵和右边源项的计算
PPE系数本质上就是压力拉普拉斯算子的离散格式中，$P_i$和$P_j$前的系数，$P_i$前的系数是主对角线元素，$P_j$前的系数是非主对角线元素。右边源项的基础是$\frac{\rho}{\Delta t}\nabla \cdot \mathbf{u}_i^*$，除此之外，压力拉普拉斯算子中壁面第二类边界条件的贡献会移项到右边成为附加源项。

为简化表达，定义

$$
L_{ij}
=
\frac{1}{r_e^2}
w_{ij}
(M_{i,3}+M_{i,4}+M_{i,5})
\mathbf{p}_{ij}
$$

以及壁面边界贡献

$$
B_i
=
\frac{1}{r_e}
\sum_{j \in Wall}
w_{ij}
(M_{i,3}+M_{i,4}+M_{i,5})
\mathbf{q}_{ij}
C_j
$$

则非自由面流体粒子的PPE装配形式为

$$
\sum_{j \in Fluid,\ j \neq i}
L_{ij}(P_j-P_i)
=
\frac{\rho}{\Delta t}
<\nabla \cdot \mathbf{u}^*>_i
- B_i
$$

对应线性系统$A\mathbf{P}=\mathbf{b}$中的元素为

$$
A_{ii}=-\sum_{j \in Fluid,\ j \neq i}L_{ij}
$$

$$
A_{ij}=L_{ij},\qquad j \in Fluid,\ j \neq i
$$

$$
b_i=
\frac{\rho}{\Delta t}
<\nabla \cdot \mathbf{u}^*>_i
- B_i
$$

如果粒子$i$不涉及壁面第二类边界条件，则$B_i=0$，且$M_i$只由流体邻居构造。

对于自由面流体粒子，直接施加压力第一类边界条件：

$$
A_{ii}=1,\qquad A_{ij}=0\ (j\neq i),\qquad b_i=P_{\mathrm{fs}}
$$

其中通常取$P_{\mathrm{fs}}=0$。因此，当前PPE求解中内部粒子、近自由面粒子以及其他非自由面流体粒子的压力拉普拉斯算子均采用type-A格式；自由面粒子不使用type-B格式，而是在压力线性系统中直接施加Dirichlet方程行。
