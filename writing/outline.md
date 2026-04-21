# 题目方向

三维空间中面向捕网反无人机场景的四智能体目标平面方形编队控制。


## 1 Introduction

### 1.1 研究背景

- 多智能体编队控制在协同拦截、协同围捕和协同捕网等任务中具有重要应用价值。
- 在捕网反无人机场景中，执行主体通常需要在三维空间中维持一个目标平面构型，以保证网口展开几何、捕获方向和包覆效果。
- 与传统二维平面编队不同，此类任务中智能体运动于三维空间，但目标构型本身是平面构型，因此不仅需要平面内形状调节，还需要处理非共面初始构型到目标平面的恢复问题。

### 1.2 相关工作

- 距离型刚性编队控制：可稳定刚性形状，但一般不直接面向目标平面恢复问题。
- 二维角度编队控制：可利用局部 bearing / angle 信息实现平面形状控制，但通常只讨论二维平面内构型。
- 3D angle rigidity 相关工作：建立了三维内角约束的刚性描述与纯方向测量控制框架，但主要面向一般三维角刚性构型，而非目标为平面构型的恢复问题。
- prescribed orientation / formation attitude 工作：控制的是整体刚性构型朝向，通常不显式处理共面恢复与无向平面法向稳定。

### 1.3 本文关注的问题

- 现有 angle-based 工作对三维目标平面编队讨论不足。
- 单纯依赖启发式角度项时，离面模态缺乏足够强的显式恢复机制。
- 当前任务中，目标法向的正反方向在几何上等价，因此不需要区分镜像和有向法向，而更关心无向平面法向意义下的目标集合恢复。

### 1.4 本文主要贡献

1. 提出一种面向三维目标平面方形编队的控制结构，在角度项和距离项之外引入共面恢复项与无向法向控制项，用于处理非共面初始构型下的目标平面恢复与平面朝向稳定。
2. 给出一个代表性离面扰动分析，说明当前启发式角度项对离面偏差仅表现出二阶敏感性，而共面误差对同一扰动是一阶敏感，从而解释显式共面恢复项的必要性。
3. 在固定分支、固定参考法向和目标构型邻域内，构造最小内部误差坐标，分析原始误差中的冗余关系，并建立闭环误差系统的局部稳定性结论。


## 2 Problem Formulation

### 2.1 Agent Dynamics

- 给出四个智能体在三维空间中的动力学模型。
- 若正文仍需兼顾实际系统，可先写理想外环模型，再说明绳索张力补偿或扰动补偿后的名义闭环分析对象。

\[
\dot p_i = v_i,\qquad
\dot v_i = u_i,\qquad i\in\{1,2,3,4\}.
\]

### 2.2 Available Measurements and Exchanged Information

- 每个智能体可在其局部坐标系下获取所需的相对 bearing 信息。
- agent 1 可测量相邻智能体的距离，并向相关智能体发送必要距离信息。
- 基于局部 bearing 与有限距离信息交换，相关智能体能够构造参考平面中心及法向控制所需量。

这里的口径统一为：

- 本文不是 purely local measurements；
- 更准确地说是 local measurements plus limited communication；
- 不要求完整全局位姿信息。

### 2.3 Desired Planar Square Formation

- 定义目标边长 \(d^\star\)。
- 定义目标平面法向 \(n_d\) 的无向等价类 \(\{n_d,-n_d\}\)。
- 定义目标平面方形构型集合。

说明：

- 本文的控制目标不是恢复唯一有向构型；
- 而是恢复目标平面方形等价类；
- 不区分由无向角度与无向法向所导致的镜像构型。

### 2.4 Assumptions

#### Assumption 1

期望编队为三维空间中的非退化平面方形构型，其目标边长为 \(d^\star\)，目标平面法向按无向量处理。

#### Assumption 2

每个智能体均可在其局部坐标系下获得控制律所需的相对 bearing 信息；agent 1 可向相邻智能体发送必要距离信息，以支持参考平面相关量的局部构造。

#### Assumption 3

每个智能体均可在其局部坐标系下获得期望法向量的表达。本文不讨论该参考量的生成与分发机制。

#### Assumption 4

系统初始构型位于目标构型附近，且用于定义角度、共面误差和法向误差的相关边向量非零、参考三角形非退化。

### 2.5 Definitions and Problem Statement

#### Definition 1: Unsigned interior angle

\[
\measuredangle ijk=\arccos\!\left(b_{ij}^\top b_{ik}\right)\in[0,\pi].
\]

该定义与三维角度刚性文献中的内角定义一致，但本文只将其作为平面内形状调节量，不直接依赖完整 3D angle rigidity 构造理论。

#### Definition 2: Coplanarity error

- 定义当前采用的共面误差，例如
\[
e_{cop}=z_{31}^\top n_{234}.
\]
- 说明该误差直接反映离面偏差。

#### Definition 3: Undirected normal error

- 定义当前平面法向；
- 定义目标法向的分支选择；
- 定义法向误差 \(e_n\)。

#### Remark 1

- 本文不讨论编队镜像区分问题。
- 原因在于：本文采用无向内角，且将目标法向视为无向量，因此 \(n_d\) 与 \(-n_d\) 代表同一目标平面朝向类。
- 因而本文研究的是目标平面方形等价类的恢复，而非唯一有向实现。

#### Problem 1

设计分布式控制律，使四智能体系统在三维空间中收敛到目标平面方形编队。具体而言，需要满足：

- 角度误差收敛至零；
- 选定距离误差收敛至零；
- 共面误差收敛至零；
- 无向法向误差收敛至零；
- 所有智能体速度收敛至零。


## 3 Main Results

### 3.1 Control Design

对每个智能体写成

\[
u_i = u_i^{ang}+u_i^{dist}+u_i^{cop}+u_i^{nor}+u_i^{vel}.
\]

各项作用如下：

- \(u_i^{ang}\)：平面内角度调节；
- \(u_i^{dist}\)：尺度与边长约束；
- \(u_i^{cop}\)：离面恢复；
- \(u_i^{nor}\)：目标平面朝向稳定；
- \(u_i^{vel}\)：速度阻尼。

#### Remark 2

- 共面项负责抑制离面模态；
- 法向项负责目标平面朝向对齐；
- 角度项与距离项负责平面内形状调节；
- 四者共同构成适用于三维目标平面编队的闭环结构。

### 3.2 Necessity of the Coplanar Recovery Term

考虑代表性离面扰动

\[
p_1=(0,d^\star,0)^T,\quad
p_2=(0,0,0)^T,\quad
p_4=(d^\star,d^\star,0)^T,\quad
p_3=(d^\star,0,h)^T,\quad |h|\ll d^\star.
\]

#### Proposition 1

在上述代表性离面扰动下：

- 三个角误差保持为零；
- 剩余角误差为 \(O(h^2)\)；
- 共面误差为 \(O(h)\)。

proof 思路：

- 逐点计算四个顶点角；
- 给出唯一被激发的角误差展开；
- 计算共面误差并作 Taylor 展开；
- 得出角通道为二阶敏感、共面通道为一阶敏感。

#### Remark 3

- Proposition 1 说明的问题，不是“最终几何上四个直角是否导向共面”，而是“角反馈通道对离面扰动是否提供足够强的一阶恢复力”。
- 这正是本文引入显式共面恢复项的原因。

### 3.3 Minimal Error Coordinates and Redundancy Analysis

- 列出原始角度误差、距离误差、共面误差、法向误差、速度误差。
- 构造最小内部误差向量。

#### Proposition 2

- 原始误差中存在冗余关系；
- 可基于独立误差构造最小内部误差坐标；
- 该最小坐标足以描述目标构型邻域内的闭环误差动力学。

#### Remark 4

- 去除冗余后，更容易识别平面内模态、离面模态和法向模态；
- 也为后续 Hurwitz 分析提供最小维度模型。

### 3.4 Local Stability Analysis

在这一部分中，单独列出局部分析条件，而不再放回 Problem Formulation：

- fixed branch；
- fixed reference normal；
- sufficiently small neighborhood；
- no switching across branch boundary。

#### Theorem 1

在给定增益条件下，最小误差系统的线性化矩阵为 Hurwitz。

proof 思路：

- 利用矩阵结构分块分析；
- 指出角度/距离模态与共面/法向模态的作用分工；
- 说明各增益项如何共同保证所有特征值位于左半平面。

#### Theorem 2

在 Assumptions 1--4 及固定分支、固定参考法向条件下，目标平面方形等价类在所考虑邻域内局部指数稳定。

proof 思路：

- Theorem 1 给出线性化系统 Hurwitz；
- 结合误差系统的光滑性和标准局部稳定性定理，推出局部指数收敛。

#### Remark 5

- 局部稳定性主要由共面项和法向项提供离面与朝向恢复；
- 角度项与距离项主要保证平面内形状；
- 该结果不直接涵盖分支切换与大范围非线性机动。

### 3.5 Simulations or Experiments

- Baseline formation recovery：从目标构型邻域内初值出发，验证整体收敛。
- Off-plane recovery experiment：采用单点抬升初值，观察共面误差、法向误差、角度误差和轨迹演化。
- Ablation study：比较完整控制律、去掉共面项、去掉法向项时的恢复性能。


## 当前建议的结果层级

- Assumption 1--4
- Definition 1--3
- Remark 1：本文不区分镜像/采用无向法向
- Problem 1
- Remark 2：控制项功能分工
- Proposition 1：代表性离面扰动下角误差与共面误差的阶次比较
- Remark 3：共面项必要性的解释
- Proposition 2：原始误差中的冗余关系
- Remark 4：最小误差坐标的作用
- Theorem 1：线性化最小误差系统的 Hurwitz 性
- Theorem 2：闭环系统局部指数稳定
- Remark 5：稳定性结论的物理解释与适用边界


## 写作上的几个统一口径

- 不把本文表述成一般 3D angle rigidity 工作；
- 将角度项称为 heuristic angle control term 或 angle-based shape regulation term；
- 将法向项表述为 undirected plane-normal control，而不是一般 rigid-body orientation control；
- 将目标集合表述为 target planar-square equivalence class，而不是唯一有向构型；
- 将信息条件表述为 local measurements plus limited communication，而不是 only local measurements。
