# 固定分支固定常值法向下四智能体无向平面编队名义最小内部误差系统的局部指数稳定性证明

## 摘要

本文研究固定分支、固定常值参考法向下，四智能体无向平面编队控制律的名义最小内部误差系统在目标方形附近的局部稳定性问题。全文首先写出当前实际外环中“名义控制 + 绳索张力 + LESO 补偿”的结构，并把绳索张力看作经 LESO 补偿的加速度通道附加项；在局部精确补偿口径下，真实外环退化为本文所分析的冻结名义外环模型。随后在代表性目标方形处构造最小内部误差坐标，证明原始误差中的一个角误差在局部上是冗余变量；再由当前控制律独立推导线性化结构矩阵，并证明相应最小内部误差系统的线性化矩阵为 Hurwitz；最后将结论推广到同一固定分支上的任意同构目标方形，并说明固定分支实现下的局部适用条件。本文不声称在任意 LESO 补偿残差、期望法向阶段切换或真正跨分支混杂轨道下都成立统一的局部指数稳定；若补偿残差不为零，则更自然的结论口径是局部 ISS、局部一致最终有界或局部 practical stability。

## 1. 问题设置

### 1.1 证明对象与实际实现对象

**定义 1（当前考虑的控制方案版本）** 为与本文证明对象保持一致，下文把“不启用饱和控制”的控制方案版本作为当前讨论对象。除饱和控制外，代码层面仍可包含 LESO、绳力支路、带迟滞记忆的无向分支选择以及阶段性法向参考生成。

与当前代码实现一致，若本文只把绳索张力并入桥接分析，则启用 LESO 时真实外环记为 $\dot p_i=v_i,\ \dot v_i=u_i^{raw}-\hat d_i+\frac{F_i^{rope}}{m}$，其中 $F_i^{rope}$ 为绳索张力合力。定义张力补偿残差 $\varepsilon_i^{rope}:=\frac{F_i^{rope}}{m}-\hat d_i$，则真实外环可等价写成 $\dot v_i=u_i^{raw}+\varepsilon_i^{rope}$。

**定义 2（本文证明对象）** 本文后续所有“局部指数稳定”结论都不直接针对定义 1 中的完整实现，而只针对上述实际外环在局部精确补偿口径 $\varepsilon_i^{rope}\equiv 0$ 下诱导出的冻结名义外环模型：控制方案中不启用饱和控制，参考法向固定为常值 $\bar n_d$，分支固定为目标分支，且同类控制项的增益统一记为常数 $k_d,\ k_\angle,\ k_v,\ k_{cop},\ k_n$；姿态内环误差等未建模附加动态不并入本文证明对象。

**假设 1（局部分析域）** 全文分析仅在目标构型附近的充分小邻域内进行，并要求所有相关边向量非零、三角形 $\{1,2,4\}$ 与 $\{2,3,4\}$ 非退化、轨道不穿越无向法向分支切换集合。

**注记 1（张力与观测补偿的桥接口径）** 若在当前工作点附近 LESO 对绳索张力的估计满足 $\hat d_i=\frac{F_i^{rope}}{m}$，则 $\varepsilon_i^{rope}\equiv 0$，真实外环恰好退化为下文的名义外环；若只知 $\varepsilon_i^{rope}$ 小、有界或渐消，则第 4 节的 Hurwitz 分析应理解为受扰系统的局部 ISS、局部一致最终有界或补偿残差渐消下的局部渐近收敛结论。

**注记 2（推导口径）** 下文中的 $J_{sq}$、$J_{raw}$、$C^\star$、$(-B^\star)$、$G^\star$、$\hat G^\star$ 及其特征多项式均从当前控制律定义和本文局部坐标独立推出，不将外部参考稿中的矩阵或谱结论作为前提。

### 1.2 冻结名义外环模型

在定义 2 下，名义外环系统写成 $\dot p_i=v_i,\ \dot v_i=u_i^{ang}+u_i^{dist}+u_i^{vel}+u_i^{cop}+u_i^{nor},\ i=1,2,3,4$。

角度项取 $u_1^{ang}=-k_\angle e_{\theta_1}(z_{12}+z_{14}),\ u_2^{ang}=-k_\angle e_{\theta_2}(z_{23}+z_{21}),\ u_3^{ang}=-k_\angle e_{\theta_3}(z_{34}+z_{32}),\ u_4^{ang}=-k_\angle e_{\theta_4}(z_{41}+z_{43})$。

距离项取 $u_1^{dist}=k_d e_{12} z_{12}+k_d e_{14} z_{14},\ u_2^{dist}=-k_d e_{12} z_{12},\ u_3^{dist}=0,\ u_4^{dist}=-k_d e_{14} z_{14}$。

速度阻尼项取 $u_{i,a}^{vel}=-k_v(v_{i,a}-v_a^\star),\ a\in\{x,y,z\}$，因此速度子块为 $-k_v I_{12}$。

共面项取 $u_3^{cop}=-k_{cop}e_{cop}n_{234},\ u_1^{cop}=u_2^{cop}=u_4^{cop}=0$，其中 $e_{cop}=z_{31}^Tn_{234}$，$n_{234}=({z_{32}\times z_{34}})/{\|z_{32}\times z_{34}\|}$。

法向项取 $u_i^{nor}=k_n(e_n\times r_i),\ i\in\{1,2,4\},\ u_3^{nor}=0$，其中 $e_n=n_{142}\times\bar n_d$，$n_{142}=({z_{14}\times z_{12}})/{\|z_{14}\times z_{12}\|}$，$r_i=p_i-(p_1+p_2+p_4)/3$。

**注记 3（共移参考系）** 文中的代表性目标方形 $p_i^\star$ 是共移坐标系或内部坐标意义下的冻结代表构型，而不是惯性系中的静止轨迹。当公共参考速度为 $v^\star$ 时，惯性系中的目标编队应理解为 $p_i(t)=p_i^\star+v^\star t+c$；因此 $\tilde v_i=v_i-v^\star$ 与各内部几何误差同时为零所描述的是匀速平移编队的稳态，而不是惯性系中的静止平衡点。

### 1.3 代表性目标方形

取固定无向平面代表分支上的代表性目标方形 $p_1^\star=[0,\ell^\star,0]^T,\ p_2^\star=[0,0,0]^T,\ p_3^\star=[\ell^\star,0,0]^T,\ p_4^\star=[\ell^\star,\ell^\star,0]^T,\ \ell^\star>0$，并取固定参考法向 $\bar n_d=[0,0,-1]^T$。

在该代表构型处有 $d_{12}^\star=d_{14}^\star=\ell^\star,\ \theta_1^\star=\theta_2^\star=\theta_3^\star=\theta_4^\star=\pi/2,\ n_{234}^\star=n_{142}^\star=\bar n_d,\ e_{cop}^\star=0,\ e_n^\star=0$，并且 $z_{12}^{\star T}z_{14}^\star=z_{21}^{\star T}z_{23}^\star=z_{32}^{\star T}z_{34}^\star=z_{43}^{\star T}z_{41}^\star=0$，$z_{14}^{\star T}z_{23}^\star=z_{12}^{\star T}z_{43}^\star=z_{21}^{\star T}z_{34}^\star=z_{32}^{\star T}z_{41}^\star=1$，从而 $z_{23}^\star=z_{14}^\star,\ z_{43}^\star=z_{12}^\star,\ z_{34}^\star=z_{21}^\star,\ z_{41}^\star=z_{32}^\star$。

## 2. 局部误差坐标与原始闭环系统

### 2.1 原始内部误差与局部构型坐标

定义原始内部误差 $\eta^{raw}=[e_{12},e_{14},e_{\theta_1},e_{\theta_2},e_{\theta_3},e_{\theta_4},e_{cop},(e_n^\perp)_1,(e_n^\perp)_2]^T\in\mathbb R^9$，其中 $e_{12}=d_{12}-\ell^\star,\ e_{14}=d_{14}-\ell^\star,\ e_{cop}=z_{31}^Tn_{234},\ e_n^\perp=E_d^T(n_{142}\times\bar n_d),\ E_d=[[1,0],[0,1],[0,0]]$，四个角误差为 $e_{\theta_1}=\angle(2,1,4)-\pi/2,\ e_{\theta_2}=\angle(3,2,1)-\pi/2,\ e_{\theta_3}=\angle(4,3,2)-\pi/2,\ e_{\theta_4}=\angle(1,4,3)-\pi/2$。

**注记 4（为什么是 8 维而不是 6 维规约）** 四个智能体在三维空间共有 12 个位置自由度。若完全不固定法向，则刚体对称通常包含 3 个平移和 3 个转动自由度；但在本文的冻结名义模型中，参考法向 $\bar n_d$ 已在惯性系中固定，因此改变 $\bar n_d$ 的两个转动自由度不再是系统对称，只剩 3 个平移自由度和 1 个绕 $\bar n_d$ 的平面内转动自由度。故内部形状自由度为 $12-4=8$。

为在目标构型附近建立局部图，取局部规约 $p_2\equiv 0$ 与 $(p_1)_x\equiv 0$，并引入 8 维局部构型坐标 $\xi=[\xi_1,\xi_2,\xi_3,\xi_4,\xi_5,\xi_6,\xi_7,\xi_8]^T$，
其中 $p_1=[0,\ell^\star+\xi_1,\xi_2]^T,\ p_3=[\ell^\star+\xi_3,\xi_4,\xi_5]^T,\ p_4=[\ell^\star+\xi_6,\ell^\star+\xi_7,\xi_8]^T$。定义原始误差映射 $h_{raw}(\xi)=\eta^{raw}$，其像记为局部可实现原始误差流形 $\mathcal M_{raw}=h_{raw}(U_\xi)\subset\mathbb R^9$。

### 2.2 原始闭环误差系统

记速度误差堆叠为 $\tilde v=[\tilde v_1^T,\tilde v_2^T,\tilde v_3^T,\tilde v_4^T]^T\in\mathbb R^{12}$，其中 $\tilde v_i=v_i-v^\star$，并定义原始闭环误差状态 $X=[(\eta^{raw})^T,\tilde v^T]^T\in\mathbb R^{21}$。

由 $\dot d_{ij}=z_{ij}^T(v_j-v_i)$、$\cos\theta_i=z_{ij}^Tz_{ik}$、$e_{cop}=z_{31}^Tn_{234}$ 与 $e_n^\perp=E_d^T(n_{142}\times\bar n_d)$ 的链式求导，定义选择矩阵 $E_i$ 使得 $\Pi_{ij}=E_j-E_i$，并记 $P_{ij}=I-z_{ij}z_{ij}^T$，$M_{ij}(p)=d_{ij}^{-1}P_{ij}\Pi_{ij}$，$S(x)$ 为叉乘矩阵，则有

$\dot\eta^{raw}=C(p)\tilde v,\ C(p)=\begin{bmatrix}C_d(p)\\ C_\theta(p)\\ C_{cop}(p)\\ C_n(p)\end{bmatrix}$，

其中

$C_d(p)=\begin{bmatrix}z_{12}^T\Pi_{12}\\ z_{14}^T\Pi_{14}\end{bmatrix}$，

$C_\theta(p)=\begin{bmatrix}-(z_{14}^TM_{12}+z_{12}^TM_{14})\\ -(z_{21}^TM_{23}+z_{23}^TM_{21})\\ -(z_{32}^TM_{34}+z_{34}^TM_{32})\\ -(z_{43}^TM_{41}+z_{41}^TM_{43})\end{bmatrix}$，

$C_{cop}(p)=n_{234}^TM_{31}+z_{31}^TN_3(p),\ N_3(p)=\frac{1}{\|z_{32}\times z_{34}\|}(I-n_{234}n_{234}^T)(-S(z_{34})M_{32}+S(z_{32})M_{34})$，

$C_n(p)=-E_d^TS(\bar n_d)N_{124}(p),\ N_{124}(p)=\frac{1}{\|z_{14}\times z_{12}\|}(I-n_{142}n_{142}^T)(-S(z_{12})M_{14}+S(z_{14})M_{12})$。

另一方面，由冻结名义控制律可得

$\dot{\tilde v}=(-B(p))\eta^{raw}-k_vI_{12}\tilde v,\ (-B(p))=[(-B_d(p)),(-B_\theta(p)),(-B_{cop}(p)),(-B_n(p))]$，

其中

$(-B_d(p))=\begin{bmatrix}k_dz_{12}&k_dz_{14}\\ -k_dz_{12}&0\\ 0&0\\ 0&-k_dz_{14}\end{bmatrix}$，

$(-B_\theta(p))=\operatorname{diag}(-k_\angle(z_{12}+z_{14}),-k_\angle(z_{23}+z_{21}),-k_\angle(z_{34}+z_{32}),-k_\angle(z_{41}+z_{43}))$，

$(-B_{cop}(p))=\begin{bmatrix}0\\0\\-k_{cop}n_{234}\\0\end{bmatrix}$，

$(-B_n(p))=\begin{bmatrix}-k_nS(r_1)E_d\\ -k_nS(r_2)E_d\\ 0_{3\times 2}\\ -k_nS(r_4)E_d\end{bmatrix}$。

因此，原始闭环误差系统可写成

$\dot X=F_{raw}(X)=\begin{bmatrix}C(p)\tilde v\\ (-B(p))\eta^{raw}-k_vI_{12}\tilde v\end{bmatrix}$，

其中 $X$ 并不在整个 $\mathbb R^{21}$ 中自由演化，而是受限于局部可实现流形 $\mathcal M_{raw}\times\mathbb R^{12}$。

在代表性目标方形处对上述原始闭环误差系统线性化，可得

$\dot X=A_{raw}^\star X=\begin{bmatrix}0_{9\times 9}&C^\star\\ (-B^\star)&-k_vI_{12}\end{bmatrix}X$

其中 $C^\star=C(p^\star)$，$(-B^\star)=(-B(p^\star))$。由第 2.2 节中的显式表达在 $p^\star$ 处代入几何关系，得到

$C^\star=\begin{bmatrix}0&1&0&0&-1&0&0&0&0&0&0&0\\ -1&0&0&0&0&0&0&0&0&1&0&0\\ \frac{1}{\ell^\star}&-\frac{1}{\ell^\star}&0&-\frac{1}{\ell^\star}&0&0&0&0&0&0&\frac{1}{\ell^\star}&0\\ -\frac{1}{\ell^\star}&0&0&\frac{1}{\ell^\star}&\frac{1}{\ell^\star}&0&0&-\frac{1}{\ell^\star}&0&0&0&0\\ 0&0&0&0&-\frac{1}{\ell^\star}&0&-\frac{1}{\ell^\star}&\frac{1}{\ell^\star}&0&\frac{1}{\ell^\star}&0&0\\ 0&\frac{1}{\ell^\star}&0&0&0&0&\frac{1}{\ell^\star}&0&0&-\frac{1}{\ell^\star}&-\frac{1}{\ell^\star}&0\\ 0&0&-\frac{\sqrt2}{2\ell^\star}&0&0&\frac{\sqrt2}{2\ell^\star}&0&0&-\frac{\sqrt2}{2\ell^\star}&0&0&\frac{\sqrt2}{2\ell^\star}\\ 0&0&-\frac{1}{\ell^\star}&0&0&\frac{1}{\ell^\star}&0&0&0&0&0&0\\ 0&0&-\frac{1}{\ell^\star}&0&0&0&0&0&0&0&0&\frac{1}{\ell^\star}\end{bmatrix}$，

$(-B^\star)=\begin{bmatrix}0&k_d&-k_\angle&0&0&0&0&0&0\\ -k_d&0&k_\angle&0&0&0&0&0&0\\ 0&0&0&0&0&0&0&\frac{k_n\ell^\star}{3}&\frac{k_n\ell^\star}{3}\\ 0&0&0&-k_\angle&0&0&0&0&0\\ k_d&0&0&-k_\angle&0&0&0&0&0\\ 0&0&0&0&0&0&0&-\frac{2k_n\ell^\star}{3}&\frac{k_n\ell^\star}{3}\\ 0&0&0&0&k_\angle&0&0&0&0\\ 0&0&0&0&-k_\angle&0&0&0&0\\ 0&0&0&0&0&0&k_{cop}&0&0\\ 0&-k_d&0&0&0&k_\angle&0&0&0\\ 0&0&0&0&0&k_\angle&0&0&0\\ 0&0&0&0&0&0&0&\frac{k_n\ell^\star}{3}&-\frac{2k_n\ell^\star}{3}\end{bmatrix}$。

进一步定义原始结构矩阵 $G^\star=C^\star(-B^\star)$，则有

$G^\star=\begin{bmatrix}-2k_d&0&k_\angle&k_\angle&0&0&0&0&0\\ 0&-2k_d&k_\angle&0&0&k_\angle&0&0&0\\ \frac{k_d}{\ell^\star}&\frac{k_d}{\ell^\star}&-\frac{2k_\angle}{\ell^\star}&\frac{k_\angle}{\ell^\star}&0&\frac{k_\angle}{\ell^\star}&0&0&0\\ \frac{k_d}{\ell^\star}&-\frac{k_d}{\ell^\star}&\frac{k_\angle}{\ell^\star}&-\frac{2k_\angle}{\ell^\star}&\frac{k_\angle}{\ell^\star}&0&0&0&0\\ -\frac{k_d}{\ell^\star}&-\frac{k_d}{\ell^\star}&0&\frac{k_\angle}{\ell^\star}&-\frac{2k_\angle}{\ell^\star}&\frac{k_\angle}{\ell^\star}&0&0&0\\ -\frac{k_d}{\ell^\star}&\frac{k_d}{\ell^\star}&\frac{k_\angle}{\ell^\star}&0&\frac{k_\angle}{\ell^\star}&-\frac{2k_\angle}{\ell^\star}&0&0&0\\ 0&0&0&0&0&0&-\frac{\sqrt2}{2\ell^\star}k_{cop}&-\frac{\sqrt2}{3}k_n&-\frac{\sqrt2}{3}k_n\\ 0&0&0&0&0&0&0&-k_n&0\\ 0&0&0&0&0&0&0&0&-k_n\end{bmatrix}$。

## 3. 冗余角误差消去与最小内部误差系统

### 3.1 冗余角误差的局部可去除性

定义最小内部误差 $\hat\eta=[e_{12},e_{14},e_{\theta_1},e_{\theta_2},e_{\theta_4},e_{cop},(e_n^\perp)_1,(e_n^\perp)_2]^T\in\mathbb R^8$，并定义最小误差映射 $\hat h(\xi)=\hat\eta$。

定义原始内部误差 $\eta^{raw}=[e_{12},e_{14},e_{\theta_1},e_{\theta_2},e_{\theta_3},e_{\theta_4},e_{cop},(e_n^\perp)_1,(e_n^\perp)_2]^T\in\mathbb R^9$，并定义原始误差映射 $h_{raw}(\xi)=\eta^{raw}$。

直接由误差定义对 $\xi$ 求导，并在 $\xi=0$ 处代入代表性目标方形，可得

$J_{sq}=D\hat h(0)=\begin{bmatrix}1&0&0&0&0&0&0&0\\ 0&0&0&0&0&1&0&0\\ -\frac{1}{\ell^\star}&0&0&0&0&0&\frac{1}{\ell^\star}&0\\ 0&0&0&-\frac{1}{\ell^\star}&0&0&0&0\\ \frac{1}{\ell^\star}&0&\frac{1}{\ell^\star}&0&0&-\frac{1}{\ell^\star}&-\frac{1}{\ell^\star}&0\\ 0&-\frac{\sqrt2}{2\ell^\star}&0&0&-\frac{\sqrt2}{2\ell^\star}&0&0&\frac{\sqrt2}{2\ell^\star}\\ 0&-\frac{1}{\ell^\star}&0&0&0&0&0&0\\ 0&-\frac{1}{\ell^\star}&0&0&0&0&0&\frac{1}{\ell^\star}\end{bmatrix}$。

由 $\det J_{sq}=-\frac{\sqrt2}{2(\ell^\star)^6}\neq 0$ 可知，存在邻域 $U_\xi$ 与 $U_{\hat\eta}$ 使得 $\hat h:U_\xi\to U_{\hat\eta}$ 为局部微分同胚，因此 $\hat\eta$ 可以作为可实现误差流形上的局部坐标。

由 $\det J_{sq}=-\frac{\sqrt2}{2(\ell^\star)^6}\neq 0$ 可知，存在邻域 $U_\xi$ 与 $U_{\hat\eta}$ 使得 $\hat h:U_\xi\to U_{\hat\eta}$ 为局部微分同胚，因此 $\hat\eta$ 也可以作为期望构型附近可实现误差集合的局部坐标。

再对原始误差映射求导可得

$J_{raw}=Dh_{raw}(0)=\begin{bmatrix}1&0&0&0&0&0&0&0\\ 0&0&0&0&0&1&0&0\\ -\frac{1}{\ell^\star}&0&0&0&0&0&\frac{1}{\ell^\star}&0\\ 0&0&0&-\frac{1}{\ell^\star}&0&0&0&0\\ 0&0&-\frac{1}{\ell^\star}&\frac{1}{\ell^\star}&0&\frac{1}{\ell^\star}&0&0\\ \frac{1}{\ell^\star}&0&\frac{1}{\ell^\star}&0&0&-\frac{1}{\ell^\star}&-\frac{1}{\ell^\star}&0\\ 0&-\frac{\sqrt2}{2\ell^\star}&0&0&-\frac{\sqrt2}{2\ell^\star}&0&0&\frac{\sqrt2}{2\ell^\star}\\ 0&-\frac{1}{\ell^\star}&0&0&0&0&0&0\\ 0&-\frac{1}{\ell^\star}&0&0&0&0&0&\frac{1}{\ell^\star}\end{bmatrix}$。

由 $J_{raw}$ 的第 5 行满足 $J_{raw}^{(5)}+J_{raw}^{(3)}+J_{raw}^{(4)}+J_{raw}^{(6)}=0$ 可知，原始 9 维误差并非独立坐标。定义图映射 $g:U_{\hat\eta}\to\mathbb R^9$ 为 $g(\hat\eta)=h_{raw}(\hat h^{-1}(\hat\eta))$，则存在唯一局部光滑函数 $\varphi:U_{\hat\eta}\to\mathbb R$ 使得

$g(\hat\eta)=[\hat\eta_1,\hat\eta_2,\hat\eta_3,\hat\eta_4,\varphi(\hat\eta),\hat\eta_5,\hat\eta_6,\hat\eta_7,\hat\eta_8]^T$。

并且其线性项满足 $e_{\theta_3}=-(e_{\theta_1}+e_{\theta_2}+e_{\theta_4})+O(\|\hat\eta\|^2)$。

因此原始闭环误差系统虽然可以写成 21 维块形式，但真正独立的内部误差只有 8 维，故后续稳定性分析应转移到最小内部误差系统上进行。

对应的一阶嵌入微分为

$S=Dg(0)=\begin{bmatrix}1&0&0&0&0&0&0&0\\ 0&1&0&0&0&0&0&0\\ 0&0&1&0&0&0&0&0\\ 0&0&0&1&0&0&0&0\\ 0&0&-1&-1&-1&0&0&0\\ 0&0&0&0&1&0&0&0\\ 0&0&0&0&0&1&0&0\\ 0&0&0&0&0&0&1&0\\ 0&0&0&0&0&0&0&1\end{bmatrix}$，

对应的选择矩阵为

$R=\begin{bmatrix}1&0&0&0&0&0&0&0&0\\ 0&1&0&0&0&0&0&0&0\\ 0&0&1&0&0&0&0&0&0\\ 0&0&0&1&0&0&0&0&0\\ 0&0&0&0&0&1&0&0&0\\ 0&0&0&0&0&0&1&0&0\\ 0&0&0&0&0&0&0&1&0\\ 0&0&0&0&0&0&0&0&1\end{bmatrix}$。

### 3.2 最小内部误差闭环系统与线性化模型

在去除冗余角误差后，定义最小闭环误差状态

$\hat X=[\hat\eta^T,\tilde v^T]^T\in\mathbb R^{20}$。

由原始闭环系统和嵌入关系 $\eta^{raw}=g(\hat\eta)$ 可得最小内部误差系统

$\dot{\hat\eta}=RC(p)\tilde v,\ \dot{\tilde v}=(-B(p))g(\hat\eta)-k_vI_{12}\tilde v$

即

$\dot{\hat X}=\hat F(\hat X)=\begin{bmatrix}RC(p(\hat\eta))\tilde v\\ (-B(p(\hat\eta)))g(\hat\eta)-k_vI_{12}\tilde v\end{bmatrix},\ \hat F(0)=0$。

在代表性目标方形处对该系统线性化，可得

$\dot{\hat X}=\hat A^\star\hat X=\begin{bmatrix}0_{8\times 8}&RC^\star\\ (-B^\star)S&-k_vI_{12}\end{bmatrix}\hat X$。

进一步定义最小结构矩阵

$\hat G^\star=RG^\star S$，

直接计算得到

$\hat G^\star=\begin{bmatrix}-2k_d&0&k_\angle&k_\angle&0&0&0&0\\ 0&-2k_d&k_\angle&0&k_\angle&0&0&0\\ \frac{k_d}{\ell^\star}&\frac{k_d}{\ell^\star}&-\frac{2k_\angle}{\ell^\star}&\frac{k_\angle}{\ell^\star}&\frac{k_\angle}{\ell^\star}&0&0&0\\ \frac{k_d}{\ell^\star}&-\frac{k_d}{\ell^\star}&0&-\frac{3k_\angle}{\ell^\star}&-\frac{k_\angle}{\ell^\star}&0&0&0\\ -\frac{k_d}{\ell^\star}&\frac{k_d}{\ell^\star}&0&-\frac{k_\angle}{\ell^\star}&-\frac{3k_\angle}{\ell^\star}&0&0&0\\ 0&0&0&0&0&-\frac{\sqrt2}{2\ell^\star}k_{cop}&-\frac{\sqrt2}{3}k_n&-\frac{\sqrt2}{3}k_n\\ 0&0&0&0&0&0&-k_n&0\\ 0&0&0&0&0&0&0&-k_n\end{bmatrix}$。

## 4. Hurwitz 判定与局部指数稳定性

### 4.1 最小线性化系统的特征方程

对第 3.2 节中的块矩阵 $\hat A^\star$ 使用分块矩阵行列式公式，可得其特征方程

$\det(\lambda I-\hat A^\star)=\det((\lambda+k_v)I_{12})\det\!\left(\lambda I_8-\frac{1}{\lambda+k_v}RC^\star(-B^\star)S\right)=(\lambda+k_v)^4\det(\lambda(\lambda+k_v)I_8-\hat G^\star)$。

另一方面，对 $\hat G^\star$ 直接计算特征多项式，可得

$\det(\lambda I-\hat G^\star)=(\lambda+k_n)^2\left(\lambda+\frac{\sqrt2}{2\ell^\star}k_{cop}\right)\left(\lambda+\frac{4k_\angle}{\ell^\star}\right)\left(\lambda^2+\frac{2(k_d\ell^\star+k_\angle)}{\ell^\star}\lambda+\frac{2k_dk_\angle}{\ell^\star}\right)^2$。

### 4.2 Hurwitz 判定

由第 4.1 节的多项式分解可知

$\sigma(\hat G^\star)=\left\{-\frac{4k_\angle}{\ell^\star},\ -\frac{k_d\ell^\star+k_\angle-\sqrt{k_d^2(\ell^\star)^2+k_\angle^2}}{\ell^\star},\ -\frac{k_d\ell^\star+k_\angle+\sqrt{k_d^2(\ell^\star)^2+k_\angle^2}}{\ell^\star},\ -\frac{k_d\ell^\star+k_\angle-\sqrt{k_d^2(\ell^\star)^2+k_\angle^2}}{\ell^\star},\ -\frac{k_d\ell^\star+k_\angle+\sqrt{k_d^2(\ell^\star)^2+k_\angle^2}}{\ell^\star},\ -\frac{\sqrt2}{2\ell^\star}k_{cop},\ -k_n,\ -k_n\right\}$。

当 $k_d>0,\ k_\angle>0,\ k_{cop}>0,\ k_n>0,\ \ell^\star>0$ 时，$\hat G^\star$ 的全部特征值严格位于开左半平面，因此 $\hat G^\star$ 为 Hurwitz 矩阵。

设 $\mu_1,\dots,\mu_8$ 为 $\hat G^\star$ 的特征值，则由第 4.1 节的特征方程有

$\det(\lambda I-\hat A^\star)=(\lambda+k_v)^4\prod_{i=1}^8(\lambda^2+k_v\lambda-\mu_i)$。

由于每个 $\mu_i<0$，故每个二次因子都可写成 $\lambda^2+k_v\lambda+|\mu_i|$，其系数全为正，因而所有根都位于开左半平面，从而 $\hat A^\star$ 为 Hurwitz 矩阵。

### 4.3 主定理

**引理 1（Hurwitz 线性化推出局部指数稳定）** 设光滑自治系统写成 $\dot x=Ax+\phi(x)$，其中 $A$ 为 Hurwitz 矩阵，且 $\|\phi(x)\|\le L\|x\|^2$。则原点局部指数稳定。

**证明** 取任意 $Q=Q^T>0$，令 $P=P^T>0$ 为 Lyapunov 方程 $A^TP+PA=-Q$ 的唯一解。定义 $V(x)=x^TPx$，则 $\lambda_{\min}(P)\|x\|^2\le V(x)\le \lambda_{\max}(P)\|x\|^2$，且 $\dot V\le -\lambda_{\min}(Q)\|x\|^2+2\|P\|L\|x\|^3$。取充分小邻域使得 $2\|P\|L\|x\|\le \lambda_{\min}(Q)/2$，即可得到 $\dot V\le -\frac{\lambda_{\min}(Q)}{2\lambda_{\max}(P)}V$，从而推出指数估计。证毕。

**定理 1（代表性目标方形处名义最小内部误差系统的局部指数稳定性）** 在定义 2、假设 1 和代表性目标方形设定下，令 $\hat X=[\hat\eta^T,\tilde v^T]^T\in\mathbb R^{20}$。若 $k_d>0,\ k_\angle>0,\ k_v>0,\ k_{cop}>0,\ k_n>0$，则最小内部误差闭环系统在原点附近可写成 $\dot{\hat X}=\hat A^\star\hat X+\hat\phi(\hat X)$，其中 $\|\hat\phi(\hat X)\|\le L\|\hat X\|^2$，且 $\hat A^\star$ 为 Hurwitz 矩阵。因此存在可实现局部邻域 $\mathcal U_{rel}$ 与常数 $M\ge 1,\ \alpha>0$，使得任意满足 $\hat X(0)\in\mathcal U_{rel}$ 的解都满足 $\|\hat X(t)\|\le Me^{-\alpha t}\|\hat X(0)\|,\ \forall t\ge 0$。

**证明** 由第 3.1 节可知，原始误差中的 $e_{\theta_3}$ 可以由 $\hat\eta$ 唯一光滑恢复，因此第 3.2 节给出的最小闭环误差系统在原点附近构成一个 20 维局部自治系统。再由第 4.2 节可知，其线性化矩阵 $\hat A^\star$ 为 Hurwitz 矩阵。于是由引理 1 可得原点局部指数稳定。证毕。

**引理 2（matched 加速度扰动下的局部 ISS / UUB 口径）** 设受扰系统写成
\[
\dot x=Ax+\phi(x)+Dw,
\]
其中 $A$ 为 Hurwitz 矩阵，$D$ 为常矩阵，且在原点附近有 $\|\phi(x)\|\le L\|x\|^2$。则存在原点附近的局部邻域与类 $\mathcal K$ 函数 $\gamma$，使得该系统关于输入 $w$ 局部输入到状态稳定。特别地，若 $w$ 有界且充分小，则状态局部一致最终有界；若 $w(t)\to 0$，则状态局部渐近收敛到原点。

**证明** 取任意 $Q=Q^T>0$，令 $P=P^T>0$ 为 Lyapunov 方程 $A^TP+PA=-Q$ 的唯一解，并定义 $V(x)=x^TPx$。则
\[
\dot V\le -\lambda_{\min}(Q)\|x\|^2+2\|P\|L\|x\|^3+2\|PD\|\|x\|\|w\|.
\]
在充分小邻域内可使 $2\|P\|L\|x\|\le \lambda_{\min}(Q)/2$，从而
\[
\dot V\le -\frac{\lambda_{\min}(Q)}{2}\|x\|^2+2\|PD\|\|x\|\|w\|.
\]
再由 Young 不等式，对任意 $a,b>0$ 有 $2ab\le \frac{\lambda_{\min}(Q)}{4}a^2+\frac{4}{\lambda_{\min}(Q)}b^2$，于是得到
\[
\dot V\le -\frac{\lambda_{\min}(Q)}{4}\|x\|^2+\frac{4\|PD\|^2}{\lambda_{\min}(Q)}\|w\|^2.
\]
再结合 $\lambda_{\min}(P)\|x\|^2\le V(x)\le \lambda_{\max}(P)\|x\|^2$ 以及标准比较引理，即得系统关于 $w$ 的局部 ISS 估计。由局部 ISS 的标准推论可知：当 $w$ 有界且充分小时，状态局部一致最终有界；当 $w(t)\to 0$ 时，状态局部渐近收敛到原点。证毕。

## 5. 坐标无关推广与固定分支实现说明

### 5.1 固定分支上任意同构目标方形的坐标无关推广

定义同一固定分支上的目标方形类 $\mathcal S(\ell^\star,\bar n_d)$ 为所有边长为 $\ell^\star$ 且平面法向代表元固定为 $\bar n_d$ 的有序正方形 $P^\dagger=(p_1^\dagger,p_2^\dagger,p_3^\dagger,p_4^\dagger)$ 的集合。对任意 $P^\dagger\in\mathcal S(\ell^\star,\bar n_d)$，存在平移向量 $b$ 与旋转矩阵 $Q\in SO(3)$ 使得 $Q\bar n_d=\bar n_d$ 且 $p_i^\dagger=Qp_i^\star+b,\ i=1,2,3,4$。

由于本文误差系统只依赖相对位置量和质心中心化量，平移向量 $b$ 在这些量中会显式相消。具体地，对任意边向量有 $q_{ij}^\dagger=p_j^\dagger-p_i^\dagger=Q(p_j^\star-p_i^\star)=Qq_{ij}^\star$，因此 $d_{ij}^\dagger=\|q_{ij}^\dagger\|=d_{ij}^\star$，$z_{ij}^\dagger=q_{ij}^\dagger/d_{ij}^\dagger=Qz_{ij}^\star$。同理，法向量满足 $n_{234}^\dagger=Qn_{234}^\star,\ n_{142}^\dagger=Qn_{142}^\star$。再注意法向控制项中的中心化向量满足 $r_i^\dagger=p_i^\dagger-(p_1^\dagger+p_2^\dagger+p_4^\dagger)/3=Q\!\left(p_i^\star-(p_1^\star+p_2^\star+p_4^\star)/3\right)=Qr_i^\star$，因此其中同样不含平移向量 $b$ 的影响。

由上述相对量和中心化量的变换关系可知，距离误差、角误差与共面误差保持不变，即 $e_{12}^\dagger=e_{12}^\star,\ e_{14}^\dagger=e_{14}^\star,\ e_{\theta_i}^\dagger=e_{\theta_i}^\star,\ e_{cop}^\dagger=e_{cop}^\star$。若切平面基取为随 $Q$ 运输的 $E_d^\dagger=QE_d$，则法向误差分量也满足 $(e_n^\perp)^\dagger=(E_d^\dagger)^T(n_{142}^\dagger\times\bar n_d)=E_d^T(n_{142}^\star\times\bar n_d)=e_n^{\perp\star}$。

同理，速度误差满足 $\tilde v^\dagger=(I_4\otimes Q)\tilde v^\star$，而冻结名义控制律中的每一项都只依赖于 $z_{ij},\ n_{234},\ n_{142},\ r_i$ 及各误差量，因此在刚体变换下满足等变性。于是有 $C^\dagger=C^\star(I_4\otimes Q^T)$，$(-B^\dagger)=(I_4\otimes Q)(-B^\star)$，从而 $G^\dagger=C^\dagger(-B^\dagger)=C^\star(-B^\star)=G^\star$，并且在最小坐标下有 $\hat G^\dagger=RG^\dagger S=\hat G^\star$，$\hat A^\dagger=T_Q\hat A^\star T_Q^{-1}$，其中 $T_Q=\operatorname{diag}(I_8,I_4\otimes Q)$。

**定理 2（固定分支上任意同构目标方形的坐标无关局部指数稳定性）** 在定理 1 的同一名义模型与增益条件下，对任意 $P^\dagger\in\mathcal S(\ell^\star,\bar n_d)$，若在 $P^\dagger$ 附近采用由 $Q$ 运输得到的局部最小内部误差坐标与切平面基，则相应最小内部误差系统的线性化矩阵 $\hat A^\dagger$ 与代表性目标方形处的 $\hat A^\star$ 共轭，从而具有相同谱。特别地，$P^\dagger$ 对应的最小内部误差平衡点也是局部指数稳定的。

**证明** 由上述刚体共轭关系可知 $\hat A^\dagger=T_Q\hat A^\star T_Q^{-1}$。由定理 1 知 $\hat A^\star$ 为 Hurwitz 矩阵，而 Hurwitz 性在相似变换下保持不变，故 $\hat A^\dagger$ 亦为 Hurwitz 矩阵。再对 $P^\dagger$ 附近的非线性系统重复定理 1 的论证，即得结论。证毕。

### 5.2 固定分支实现的局部说明

本文直接假设证明对象中不启用饱和控制，因此名义模型中的 $k_d,\ k_\angle,\ k_v,\ k_{cop},\ k_n$ 从一开始就是常增益参数，而不是由饱和函数在小邻域内退化得到。与此同时，当前实际实现在启用 LESO 时满足 $\dot p_i=v_i,\ \dot v_i=u_i^{raw}-\hat d_i+\frac{F_i^{rope}}{m}=u_i^{raw}+\varepsilon_i^{rope}$，其中 $\varepsilon_i^{rope}:=\frac{F_i^{rope}}{m}-\hat d_i$。故本小节需要说明的不是“把绳力排除在外”，而是“在当前工作点附近，绳索张力被 LESO 补偿后，真实外环如何退化为本文的固定分支冻结名义模型”。

**命题 1（LESO 对绳索张力的局部补偿所诱导的名义退化）** 若在某局部邻域内对所有 $i\in\{1,2,3,4\}$ 都有 $\hat d_i=\frac{F_i^{rope}}{m}$，则真实外环动力学与第 1.2 节给出的冻结名义外环完全一致。

**证明** 由定义 $\varepsilon_i^{rope}=\frac{F_i^{rope}}{m}-\hat d_i$ 可知，在该局部邻域内 $\varepsilon_i^{rope}\equiv 0$。于是实际外环速度方程化为 $\dot v_i=u_i^{raw}$，这与第 1.2 节的冻结名义外环逐项一致，故结论成立。证毕。

实际实现中的无向法向参考不是固定常值，而是通过迟滞分支逻辑选择 $\bar n_d=s_kn_d$，其中 $s_k\in\{+1,-1\}$ 由对齐量 $c_k=n_{142}(t_k)^Tn_d$ 与阈值 $\delta_h\in(0,1)$ 更新，即 $s_k=+1$ 当 $c_k>\delta_h$，$s_k=-1$ 当 $c_k<-\delta_h$，其余情形保持上一时刻分支。

**命题 2（目标分支附近的局部无切换性）** 设 $P^\dagger\in\mathcal S(\ell^\star,\bar n_d)$，并令与该固定分支对应的符号为 $s^\star\in\{+1,-1\}$，满足 $\bar n_d=s^\star n_d$。若初始分支已经取为正确值 $s_0=s^\star$，则存在 $P^\dagger$ 附近的局部邻域 $\mathcal U_{sw}$，使得只要轨道停留在 $\mathcal U_{sw}$ 内，迟滞分支状态始终保持 $s_k=s^\star$，从而实际分支逻辑在该邻域内退化为固定分支逻辑。

**证明** 在 $P^\dagger$ 处有 $n_{142}^\dagger=\bar n_d=s^\star n_d$，故 $s^\star n_{142}^{\dagger T}n_d=1$。由 $n_{142}(p)$ 的连续性，存在平衡点邻域 $\mathcal U_{sw}$ 使得 $|s^\star n_{142}(p)^Tn_d-1|<\varepsilon$，其中 $0<\varepsilon<1-\delta_h$。于是对任意 $p\in\mathcal U_{sw}$ 都有 $s^\star n_{142}(p)^Tn_d>1-\varepsilon>\delta_h$。由于初始分支已满足 $s_0=s^\star$，迟滞更新判据在全部更新时刻都会保持在目标分支一侧，因此不会触发跨分支切换，故 $s_k$ 始终保持为 $s^\star$。证毕。

在命题 2 给出的固定分支局部邻域内，即使 LESO 对绳索张力的补偿不精确，真实最小内部误差系统仍可写成受扰形式
\[
\dot{\hat X}=\hat A^\star \hat X+\hat\phi(\hat X)+D_\varepsilon \varepsilon^{rope},
\qquad
D_\varepsilon=
\begin{bmatrix}
0_{8\times 12}\\
I_{12}
\end{bmatrix},
\]
其中 $\varepsilon^{rope}=[(\varepsilon_1^{rope})^T,(\varepsilon_2^{rope})^T,(\varepsilon_3^{rope})^T,(\varepsilon_4^{rope})^T]^T\in\mathbb R^{12}$ 为 LESO 与绳索张力之间的补偿残差堆叠向量，且 $\hat\phi(\hat X)$ 满足与定理 1 相同的局部二阶有界性质。

**命题 3（LESO 与绳索张力存在补偿残差时的局部 ISS / UUB 口径）** 在定理 1 与命题 2 的条件下，若 $\varepsilon^{rope}$ 在某局部邻域内有界，则上述真实最小内部误差系统关于输入 $\varepsilon^{rope}$ 局部输入到状态稳定。特别地：

1. 若 $\|\varepsilon^{rope}(t)\|$ 一致有界且充分小，则 $\hat X(t)$ 局部一致最终有界；
2. 若 $\varepsilon^{rope}(t)\to 0$，则 $\hat X(t)\to 0$；
3. 若进一步有 $\varepsilon^{rope}\equiv 0$，则恢复定理 1 和定理 2 的局部指数稳定结论。

**证明** 由命题 2 可知，在该局部邻域内实际分支逻辑退化为固定分支逻辑；因此补偿残差仅通过速度通道进入最小内部误差系统，并可写成上述受扰形式。又由于定理 1 已证明 $\hat A^\star$ 为 Hurwitz，且 $\hat\phi(\hat X)$ 满足局部二阶有界条件，故直接对受扰系统应用引理 2，即得该系统关于 $\varepsilon^{rope}$ 的局部 ISS 结论。由局部 ISS 的标准推论分别得到局部一致最终有界、残差渐消下的局部渐近收敛以及零残差下恢复局部指数稳定。证毕。

由命题 1 与命题 2 可知，在“控制方案中不启用饱和控制、LESO 在该局部邻域内精确补偿绳索张力、初始分支已正确、轨道始终停留在局部邻域内”这一组条件下，实际外环实现退化为本文的固定分支冻结名义模型，因此定理 2 与推论 1 可直接作用于该实现。若 LESO 与绳索张力之间存在补偿残差，即 $\varepsilon_i^{rope}\not\equiv 0$，则不能再直接声称严格的局部指数稳定，此时更自然的结论口径如命题 3 所示，应理解为局部 ISS、局部一致最终有界或残差渐消下的局部渐近收敛。

## 6. 物理可解释误差的回传结论

**推论 1（角误差、距离误差、速度误差的局部指数收敛）** 在定理 2 的条件下，除最小误差状态 $\hat X=[\hat\eta^T,\tilde v^T]^T$ 外，常用的物理误差量也都局部指数收敛。具体而言，存在 $P^\dagger$ 附近的可实现局部邻域与常数 $C>0,\ \alpha>0$，使得四个角误差 $e_{\theta_1},e_{\theta_2},e_{\theta_3},e_{\theta_4}$、独立边长误差 $e_{12}=d_{12}-\ell^\star,\ e_{14}=d_{14}-\ell^\star,\ e_{23}=d_{23}-\ell^\star,\ e_{34}=d_{34}-\ell^\star$、两条对角线误差 $e_{13}=d_{13}-\sqrt2\ell^\star,\ e_{24}=d_{24}-\sqrt2\ell^\star$、共面误差 $e_{cop}$、法向误差分量 $e_n^\perp$ 与速度误差 $\tilde v_i-v^\star,\ i=1,2,3,4$ 均满足 $|e(t)|\le Ce^{-\alpha t}\|\hat X(0)\|$。

**证明** 在 $P^\dagger$ 附近，$\hat h$ 是局部微分同胚，因此局部构型坐标 $\xi=\hat h^{-1}(\hat\eta)$ 是 $\hat\eta$ 的光滑函数。四个角误差、所有独立边长误差、对角线误差、$e_{cop}$ 与 $e_n^\perp$ 都是 $\xi$ 的光滑函数，并在目标方形处取零值；其中 $e_{\theta_3}$ 还满足 $e_{\theta_3}=\varphi(\hat\eta)$。由均值定理，在充分小邻域内对每个上述误差 $e$ 都存在常数 $c_e>0$ 使得 $|e|\le c_e\|\hat\eta\|$。再结合定理 2 的估计 $\|\hat\eta(t)\|+\|\tilde v(t)\|\le Me^{-\alpha t}\|\hat X(0)\|$，即可得出所有这些物理误差都以同一指数阶收敛到零。证毕。

## 7. 结论与适用边界

本文在固定无向平面代表分支、固定常值参考法向、且控制方案中不启用饱和控制的冻结名义模型口径下，构造了四智能体方形编队的最小内部误差坐标，并通过 $J_{raw}$ 的显式行关系证明原始误差中的 $e_{\theta_3}$ 在局部上是冗余变量。随后，本文由当前控制律独立推导了 $C^\star$、$(-B^\star)$、$G^\star$ 与 $\hat G^\star$，证明了最小内部误差系统的线性化矩阵 $\hat A^\star$ 为 Hurwitz 矩阵，从而建立了代表性目标方形处名义最小内部误差系统的局部指数稳定性。进一步地，借助刚体共轭关系，该结论被推广到同一固定分支上的任意同构目标方形；若实际实现中采用与本文一致的无饱和固定分支控制方案，并且 LESO 在相应局部邻域内将绳索张力准确补偿，则真实外环退化为相应冻结外环，因而上述结论可直接作用于该实现。

需要强调的是，本文当前结论支持的最强口径是：“固定分支、固定常值法向、无饱和控制、且 LESO 对绳索张力实现局部精确补偿”条件下的最小内部误差系统局部指数稳定。若 LESO 张力补偿残差不为零，则更合适的结论应为局部 ISS、局部一致最终有界或局部 practical stability，而不是统一的严格局部指数稳定。重新引入饱和控制、姿态内环误差、期望法向阶段切换或连续变化，以及真正跨越切换带的混杂稳定性问题，仍属于后续工作。

## 附录 A. 关键矩阵条目的复核路径

为便于审稿和答辩复核，本文关键矩阵的中间求导路径归纳如下。

1. 距离误差条目来自 $d_{ij}=\|p_j-p_i\|$ 的一阶微分，即 $D d_{ij}=z_{ij}^T\Pi_{ij}$。

2. 角误差条目来自 $\cos\theta_i=z_{ij}^Tz_{ik}$ 的链式求导，即 $D e_{\theta_i}=-(z_{ik}^TM_{ij}+z_{ij}^TM_{ik})$，并在代表性目标方形处使用 $\sin\theta_i^\star=1$。

3. 共面误差条目来自 $e_{cop}=z_{31}^Tn_{234}$ 的链式求导，即 $D e_{cop}=n_{234}^TM_{31}+z_{31}^TN_3(p)$。

4. 法向误差条目来自 $e_n^\perp=E_d^T(n_{142}\times\bar n_d)$ 的链式求导，即 $D e_n^\perp=-E_d^TS(\bar n_d)N_{124}(p)$。

5. $J_{sq}$ 与 $J_{raw}$ 由第 2 节中定义的映射 $\hat h(\xi)$ 与 $h_{raw}(\xi)$ 对局部坐标 $\xi$ 直接求导得到；$C^\star$ 则由本附录前四条在 $p^\star$ 处代入几何条件得到；$(-B^\star)$ 由冻结名义控制律逐项在 $p^\star$ 处代入得到；$G^\star$ 与 $\hat G^\star$ 分别由矩阵乘积 $C^\star(-B^\star)$ 和 $RG^\star S$ 给出。

6. 配套的符号核验脚本位于 `analysis/derive_minimal_linearization.py`，可用于独立复核 $J_{sq}$、$J_{raw}$、$C^\star$、$(-B^\star)$、$G^\star$、$\hat G^\star$ 与特征多项式。




# 分支选择

## 1. 基础变量与候选集合
设定当前实际的法向量为 **\(n\)** ，期望的（未确定符号的）基准方向为 **\(n_d\)**。
首先构建候选集合，选择与当前 \(n\) 误差最小的方向：
\[
 \{n_d, -n_d\} 
\]

## 2. 计算点积判断方向差距
通过点积来衡量当前实际方向 \(n\) 与期望基准方向 \(n_d\) 之间的对齐程度：
\[
 c = n^\top n_d 
\]
> *注：\(c\) 的大小接近 1 表示两向量非常接近，接近 -1 表示反向，接近 0 表示接近垂直。*

## 3. 面临的难题
如果单纯通过判断 \(c\) 的正负来确定最终的 \(n_d\) 方向（即 \(c>0\) 选 \(n_d\)，\(c<0\) 选 \(-n_d\)），则会遇到以下问题：
> **问题：** 如果编队要求垂直方向排列，此时 \(c\) 会退化为 0（或接近0）。在这种情况下，极小的噪声就会导致正负号频繁跳变，从而引发**控制器剧烈振荡和不稳定**。

## 4. 改进的控制方法（带记忆/迟滞的符号切换逻辑）
为了避免上述问题，引入一个状态变量 \(s_k\)（符号选择因子），并设置一个死区阈值 \(\delta\)。核心的符号选择逻辑如下：

$$
s_k = 
\begin{cases} 
+1, & c > \delta \\
-1, & c < -\delta \\
s_{k-1}, & |c| \le \delta
\end{cases}
$$

**逻辑解析：**
*   **当 \(c > \delta\) 时**：当前向量与 \(n_d\) 方向偏差很小，选择 \(+n_d\)。
*   **当 \(c < -\delta\) 时**：当前向量几乎与 \(n_d\) 反向，选择 \(-n_d\)。
*   **当 \(|c| \le \delta\) 时**（即处于“垂直”或“模糊死区”）：**保持上一时刻的选择方向** (\(s_{k-1}\))，不做改变。

## 5. 最终确定的方向
根据计算出的 \(s_k\)，得到用于控制器的最终期望法向量 \(\bar{n}_d\)：
\[
 \bar{n}_d = s_k n_d 
\]

---

## 💡 方法总结
这种控制方法的核心思想是**将“定向控制”转化为“无向控制”**。通过设置死区阈值 \(\delta\) 并结合上一次状态 \(s_{k-1}\) 来形成迟滞（Hysteresis）特性。只有当点积 \(c\) 的绝对值足够大、明确偏离垂直方向时，才改变符号；在垂直方向附近，强制保持状态。这有效消除了控制器在编队垂直姿态附近的抖动和啸叫现象。