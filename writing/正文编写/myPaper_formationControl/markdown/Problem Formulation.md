在这篇文章中我们解决的是四智能体编队曳引弹性网在3D空间中完成反无人机过程中保持编队刚性的问题。

将四个智能体视为在三维空间中运动的双积分系统。考虑系统不确定性以及外部时变扰动，第 \(i\) 个智能体的动力学可写为

\[
\left\{
\begin{aligned}
\dot p_i(t)&=v_i(t),\\
\dot v_i(t)&=u_i(t)+f_i(t)+d_i(t),
\end{aligned}
\right.
\qquad i=1,2,3,4.
\tag{1}
\]

其中，\(p_i(t)\in\mathbb R^3\) 和 \(v_i(t)\in\mathbb R^3\) 分别表示第 \(i\) 个智能体的位置和速度，\(u_i(t)\in\mathbb R^3\) 为控制输入，\(f_i(t)\in\mathbb R^3\) 表示系统不确定性和未建模动力学，\(d_i(t)\in\mathbb R^3\) 表示外部时变扰动。在反无人机过程中，弹性网的张力对agent个体的扰动是不可忽略的，因此本文将agent受到弹性网的张力并入 \(d_i(t)\) 中统一处理。将上述不确定项合并为集总扰动

\[
D_i(t)=f_i(t)+d_i(t).
\tag{2}
\]

于是系统动力学可等价写为

\[
\left\{
\begin{aligned}
\dot p_i(t)&=v_i(t),\\
\dot v_i(t)&=u_i(t)+D_i(t),
\end{aligned}
\right.
\qquad i=1,2,3,4.
\tag{3}
\]

为了后续扰动补偿，记由扰动观测器给出的集总扰动估计值为 \(\hat D_i(t)\)，并定义估计误差为

\[
\tilde D_i(t)=\hat D_i(t)-D_i(t).
\]

具体观测器结构将在控制律设计部分给出。

对每一个无人机都做如下假设:
##### Assumption1 
每一个agent都能通过LOS获取控制律所需的相对方向信息
##### Assumption2(存疑是否要加入)
每个agent都可以获得期望公共速度还有期望法向量方向在其局部坐标系下的表示。
##### Assumption3
编队初始构型在期望构型的领域内，每一个智能体的位置都是独立的并且非共线。
##### Assumption4
集总扰动 \(D_i(t)=f_i(t)+d_i(t)\) 有界，其中 \(d_i(t)\) 包含外部环境扰动以及捕网绳索张力对第 \(i\) 个 UAV 产生的等效影响。


本文考虑的是四智能体捕网编队，顶点集合可以描述为

\[
\mathcal V=\{1,2,3,4\}.
\]

其中，第 \(i\) 个智能体的位置为 \(p_i\in\mathbb R^3\)，整体构型记为

\[
p=\left[p_1^T,p_2^T,p_3^T,p_4^T\right]^T\in\mathbb R^{12}.
\]

对任意 \(i,j\in\mathcal V\)，定义相对位置、相对距离和 bearing 向量为

\[
q_{ij}=p_j-p_i,\qquad
d_{ij}=\|q_{ij}\|,\qquad
b_{ij}=\frac{q_{ij}}{\|q_{ij}\|}.
\]

在Assumption 1下，agent #i 可以获取 \(b_{ij}\)、\(b_{ik}\) 的方向。因此，可以进一步由局部 bearing 信息构造相应的内角大小 \(\alpha_i\)，其中 \(\alpha_{jik}\in[0,\pi]\)。\(\alpha_{jik}\) 可以表示为：

\[
\alpha_{jik}=\arccos\left(b_{ij}^{T}b_{ik}\right).
\]

角度集合表示捕网的平面内形状约束为：

\[
\mathcal A=\{(2,1,4),(3,2,1),(4,3,2),(1,4,3)\}.
\]

其中，每个三元组 \((j,i,k)\in\mathcal A\) 表示以 UAV \(i\) 为顶点、由 UAV \(j\)、UAV \(i\)、UAV \(k\) 形成的内角。

对于期望平面正方形捕网，角度约束要求

\[
\alpha_{jik}^*(t)=\frac{\pi}{2},
\qquad
(j,i,k)\in\mathcal A.
\]

角度约束可以刻画期望形状，但由于角度在整体平移、旋转和尺度变化下保持不变，仅依靠角度约束不能确定编队的实际尺寸(scale)。
在本文中采用角度-距离结合的几何控制方法来刻画描述捕网的平面构型。针对平面正方形捕网这种四个智能体共圆的情况，引入距离约束边集ref[Multi-agent distributed formation control using local information measurements]。

\[
\mathcal E_d=\{(1,2),(1,4)\}.
\]

结合直角约束，这两个距离约束能够在局部确定捕网正方形的尺度。对应的尺度约束为

\[
d_{ij}=d_{ij}^*(t),
\qquad
(i,j)\in\mathcal E_d,
\]

其中 \(d_{ij}^*(t)=\ell^*\)，\(\ell^*>0\) 为期望捕网边长，\(d_{ij}=\|p_j-p_i\|\)。

<!-- 进一步地，本文考虑的是三维空间中的平面捕网编队。角度约束和距离约束只能刻画平面内的局部形状与尺度，不能单独保证四个 UAV 在机动过程中始终恢复到同一平面内。因此，需要额外引入共面约束。至于共面项的重要性说明将在下一节给出。 -->

进一步地，本文考虑三维空间中的平面捕网编队。需要考虑共面保持的问题，因此需要额外引入共面约束项。至于共面项的重要性将在下一节说明



\[
n_{234}=\frac{b_{32}\times b_{34}}
{\|b_{32}\times b_{34}\|},
\qquad
e_{\mathrm{cop}}=b_{31}^{T}n_{234}.
\]

fig1

如fig1所示，当 \(e_{\mathrm{cop}}=0\) 时 agent# 1 位于由 agent# 2、agent# 3 和 agent# 4 确定的平面内，从而四个 agent 共面。

为了增加捕获无人机的成功率，在本文中我们选择显式控制编队的法向量，我们选择agent# 1 agent# 2 agent# 4构成的平面为参考平面，其法向量为：

\[
n_{142}=\frac{b_{14}\times b_{12}}
{\|b_{14}\times b_{12}\|}.
\]

在捕网任务中，网面的正面和反面具有物理等价性，因此 \(n_d\) 和 \(-n_d\) 表示同一个期望捕获平面。本文采用无向法向量目标集合

\[
\mathcal N_d=\{n_d,-n_d\},
\]

并定义期望法向量为代表元$\bar n_d$


因此，定义法向量误差为：

\[
e_n = n_{142} \times \bar n_d
\]

当 \(e_n=0\) 时，有 \(n_{142}=n_d\) 或 \(n_{142}=-n_d\)，即当前捕网平面与期望捕获平面在无向意义下一致。

综上，本文期望捕网构型可由角度集合、距离边集、共面约束和无向法向量约束共同描述为

针对扰动和内部误差

\[
\mathcal S(\ell^*,n_d)=\left\{
p\in\mathbb R^{12}
\ \middle|\
\begin{aligned}
&\alpha_{jik}=\alpha_{jik}^*(t),\quad (j,i,k)\in\mathcal A,\\
&d_{ij}=d_{ij}^*(t),\quad (i,j)\in\mathcal E_d,\\
&e_{\mathrm{cop}}=0,\\
&n_{142}=\bar n_d
\end{aligned}
\right\}.
\]

其中，角度集合 \(\mathcal A\) 用于保持捕网的正方形形状，距离边集 \(\mathcal E_d\) 用于确定捕网尺度，共面约束用于保证四个顶点形成真实的物理网面，无向法向量约束用于描述捕网平面的空间朝向。

<!-- 对于一般三维非共面刚体编队，reflection ambiguity 会产生关于参考平面对称的两个不同空间构型；但对于本文考虑的共面捕网方形编队，关于编队平面本身的反射不会改变四个 UAV 的位置。因此，本文中的“正反侧歧义”主要表现为平面法向量 n 与 −n 的等价性，而不是两个不同点构型之间的歧义。 -->

为了获取必要的距离信息，做出下方假设
##### Assumption5
有且仅有Agent #1可以测量与其相邻的两个智能体的距离，并且通信给相邻智能体。


##### Problem 1

设计一种适用于四 UAV 捕网任务的编队控制律及扰动观测补偿机制，使得在存在系统不确定性、外部扰动以及捕网绳索张力等效影响的情况下，四个智能体能够在三维空间中形成并保持期望平面正方形捕网构型，同时实现捕网平面法向量控制和公共速度跟踪。该 formation maneuvering 任务由三个子目标组成：a) planar square shape control；b) plane-normal control；c) velocity tracking。

要求闭环系统从期望构型附近出发时满足如下收敛目标。

\[
\lim_{t\to+\infty}
\left(\alpha_{jik}(t)-\alpha_{jik}^*(t)\right)=0,
\qquad
(j,i,k)\in\mathcal A.
\]

\[
\lim_{t\to+\infty}
\left(d_{ij}(t)-d_{ij}^*(t)\right)=0,
\qquad
(i,j)\in\mathcal E_d.
\]


\[
\lim_{t\to+\infty}e_{\mathrm{cop}}(t)=0.
\]


\[
\lim_{t\to+\infty}
\left(e_{n}(t)\right)=0.
\]


\[
\lim_{t\to+\infty}
\left(v_i(t)-v_c^*(t)\right)=0,
\qquad
i\in\mathcal V.
\]

<!-- 当观测器估计误差 \(\tilde D_i(t)\) 有界时，要求上述几何误差和速度跟踪误差保持有界；若 \(\tilde D_i(t)\to0\)，则闭环系统恢复上述理想补偿情形下的编队机动目标。 -->