# 四智能体无向平面编队二阶积分最小内部误差系统的局部指数稳定性证明

## 摘要

本文研究四智能体无向平面编队二阶积分闭环模型在目标方形附近的局部指数稳定性问题。证明对象包含二阶积分器动力学、距离项、角度项、速度阻尼项、共面项和法向项。全文先给出闭环模型，再在代表性目标方形处构造最小内部误差坐标，证明原始误差中的一个角误差在局部上是冗余变量；随后由控制律直接推导线性化结构矩阵，并证明相应最小内部误差系统的线性化矩阵为 Hurwitz；最后说明代表性目标方形处的结论可推广到任意同构目标方形。

## 1. 问题设置

### 1.1 证明对象

**定义 1（二阶积分闭环模型）** 本文只考虑四智能体的二阶积分器动力学

$\dot p_i=v_i,\qquad \dot v_i=u_i,\qquad i=1,2,3,4,$

其中 $p_i\in\mathbb R^3$、$v_i\in\mathbb R^3$，控制输入 $u_i$ 由角度项、距离项、速度阻尼项、共面项和法向项组成。参考法向 $\bar n_d$ 在局部线性化中作为给定常量处理，所有控制增益 $k_d,k_\angle,k_v,k_{cop},k_n$ 为正常数，且不启用输入饱和。

**假设 1（局部分析域）** 全文分析仅在目标构型附近的充分小邻域内进行，并要求所有相关边向量非零，三角形 $\{1,2,4\}$ 与 $\{2,3,4\}$ 非退化，从而所用误差映射与闭环向量场在该邻域内光滑。

**注记 1（推导口径）** 下文中的 $J_{sq}$、$J_{raw}$、$C^\star$、$(-B^\star)$、$G^\star$、$\hat G^\star$ 及其特征多项式均从本文二阶积分闭环控制律和局部坐标直接推出，不把外部参考稿中的矩阵或谱结论作为前提。

### 1.2 二阶积分闭环模型

在定义 1 下，闭环二阶积分模型写成 $\dot p_i=v_i,\ \dot v_i=u_i^{ang}+u_i^{dist}+u_i^{vel}+u_i^{cop}+u_i^{nor},\ i=1,2,3,4$。

角度项取 $u_1^{ang}=-k_\angle e_{\theta_1}(z_{12}+z_{14}),\ u_2^{ang}=-k_\angle e_{\theta_2}(z_{23}+z_{21}),\ u_3^{ang}=-k_\angle e_{\theta_3}(z_{34}+z_{32}),\ u_4^{ang}=-k_\angle e_{\theta_4}(z_{41}+z_{43})$。

距离项取 $u_1^{dist}=k_d e_{12} z_{12}+k_d e_{14} z_{14},\ u_2^{dist}=-k_d e_{12} z_{12},\ u_3^{dist}=0,\ u_4^{dist}=-k_d e_{14} z_{14}$。

速度阻尼项取 $u_{i,a}^{vel}=-k_v(v_{i,a}-v_a^\star),\ a\in\{x,y,z\}$，因此速度子块为 $-k_v I_{12}$。

共面项取 $u_3^{cop}=-k_{cop}e_{cop}n_{234},\ u_1^{cop}=u_2^{cop}=u_4^{cop}=0$，其中 $e_{cop}=z_{31}^Tn_{234}$，$n_{234}=({z_{32}\times z_{34}})/{\|z_{32}\times z_{34}\|}$。

法向项取 $u_i^{nor}=k_n(e_n\times r_i),\ i\in\{1,2,4\},\ u_3^{nor}=0$，其中 $e_n=n_{142}\times\bar n_d$，$n_{142}=({z_{14}\times z_{12}})/{\|z_{14}\times z_{12}\|}$，$r_i=p_i-(p_1+p_2+p_4)/3$。

**注记 2（共移参考系）** 下文中的代表性目标方形 $p_i^\star$ 是共移坐标系或内部坐标意义下的冻结代表构型，而不是惯性系中的静止轨迹。当公共参考速度为 $v^\star$ 时，惯性系中的目标编队应理解为 $p_i(t)=p_i^\star+v^\star t+c$；因此 $\tilde v_i=v_i-v^\star$ 与各内部几何误差同时为零所描述的是匀速平移编队的稳态，而不是惯性系中的静止平衡点。

### 1.3 代表性目标方形

取如下代表性目标方形 $p_1^\star=[0,\ell^\star,0]^T,\ p_2^\star=[0,0,0]^T,\ p_3^\star=[\ell^\star,0,0]^T,\ p_4^\star=[\ell^\star,\ell^\star,0]^T,\ \ell^\star>0$，并取 $\bar n_d=[0,0,-1]^T$。

在该代表构型处有 $d_{12}^\star=d_{14}^\star=\ell^\star,\ \theta_1^\star=\theta_2^\star=\theta_3^\star=\theta_4^\star=\pi/2,\ n_{234}^\star=n_{142}^\star=\bar n_d,\ e_{cop}^\star=0,\ e_n^\star=0$，并且 $z_{12}^{\star T}z_{14}^\star=z_{21}^{\star T}z_{23}^\star=z_{32}^{\star T}z_{34}^\star=z_{43}^{\star T}z_{41}^\star=0$，$z_{14}^{\star T}z_{23}^\star=z_{12}^{\star T}z_{43}^\star=z_{21}^{\star T}z_{34}^\star=z_{32}^{\star T}z_{41}^\star=1$，从而 $z_{23}^\star=z_{14}^\star,\ z_{43}^\star=z_{12}^\star,\ z_{34}^\star=z_{21}^\star,\ z_{41}^\star=z_{32}^\star$。


## 2. 局部误差坐标与原始闭环系统

### 2.1 原始内部误差与局部构型坐标

定义原始内部误差 $\eta^{raw}=[e_{12},e_{14},e_{\theta_1},e_{\theta_2},e_{\theta_3},e_{\theta_4},e_{cop},(e_n^\perp)_1,(e_n^\perp)_2]^T\in\mathbb R^9$，其中 $e_{12}=d_{12}-\ell^\star,\ e_{14}=d_{14}-\ell^\star,\ e_{cop}=z_{31}^Tn_{234},\ e_n^\perp=E_d^T(n_{142}\times\bar n_d),\ E_d=[[1,0],[0,1],[0,0]]$，四个角误差为 $e_{\theta_1}=\angle(2,1,4)-\pi/2,\ e_{\theta_2}=\angle(3,2,1)-\pi/2,\ e_{\theta_3}=\angle(4,3,2)-\pi/2,\ e_{\theta_4}=\angle(1,4,3)-\pi/2$。

**注记 3（局部坐标维数）** 在给定 $\bar n_d$ 的局部模型中，平移自由度和绕 $\bar n_d$ 的平面内转动自由度不影响内部误差。故可用 8 维局部构型坐标描述目标构型附近的独立内部几何变化。

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

另一方面，由上述二阶积分闭环控制律可得

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

**定理 1（代表性目标方形处二阶积分最小内部误差系统的局部指数稳定性）** 在定义 1、假设 1 和代表性目标方形设定下，令 $\hat X=[\hat\eta^T,\tilde v^T]^T\in\mathbb R^{20}$。若 $k_d>0,\ k_\angle>0,\ k_v>0,\ k_{cop}>0,\ k_n>0$，则最小内部误差闭环系统在原点附近可写成 $\dot{\hat X}=\hat A^\star\hat X+\hat\phi(\hat X)$，其中 $\|\hat\phi(\hat X)\|\le L\|\hat X\|^2$，且 $\hat A^\star$ 为 Hurwitz 矩阵。因此存在可实现局部邻域 $\mathcal U_{rel}$ 与常数 $M\ge 1,\ \alpha>0$，使得任意满足 $\hat X(0)\in\mathcal U_{rel}$ 的解都满足 $\|\hat X(t)\|\le Me^{-\alpha t}\|\hat X(0)\|,\ \forall t\ge 0$。

**证明** 由第 3.1 节可知，原始误差中的 $e_{\theta_3}$ 可以由 $\hat\eta$ 唯一光滑恢复，因此第 3.2 节给出的最小闭环误差系统在原点附近构成一个 20 维局部自治系统。再由第 4.2 节可知，其线性化矩阵 $\hat A^\star$ 为 Hurwitz 矩阵。于是由引理 1 可得原点局部指数稳定。证毕。

## 5. 坐标无关推广

### 5.1 任意同构目标方形的坐标无关推广

定义目标方形类 $\mathcal S(\ell^\star,\bar n_d)$ 为所有边长为 $\ell^\star$ 且平面法向为 $\bar n_d$ 的有序正方形 $P^\dagger=(p_1^\dagger,p_2^\dagger,p_3^\dagger,p_4^\dagger)$ 的集合。对任意 $P^\dagger\in\mathcal S(\ell^\star,\bar n_d)$，存在平移向量 $b$ 与旋转矩阵 $Q\in SO(3)$ 使得 $Q\bar n_d=\bar n_d$ 且 $p_i^\dagger=Qp_i^\star+b,\ i=1,2,3,4$。

由于本文误差系统只依赖相对位置量和质心中心化量，平移向量 $b$ 在这些量中会显式相消。具体地，对任意边向量有 $q_{ij}^\dagger=p_j^\dagger-p_i^\dagger=Q(p_j^\star-p_i^\star)=Qq_{ij}^\star$，因此 $d_{ij}^\dagger=\|q_{ij}^\dagger\|=d_{ij}^\star$，$z_{ij}^\dagger=q_{ij}^\dagger/d_{ij}^\dagger=Qz_{ij}^\star$。同理，法向量满足 $n_{234}^\dagger=Qn_{234}^\star,\ n_{142}^\dagger=Qn_{142}^\star$。再注意法向控制项中的中心化向量满足 $r_i^\dagger=p_i^\dagger-(p_1^\dagger+p_2^\dagger+p_4^\dagger)/3=Q\!\left(p_i^\star-(p_1^\star+p_2^\star+p_4^\star)/3\right)=Qr_i^\star$，因此其中同样不含平移向量 $b$ 的影响。

由上述相对量和中心化量的变换关系可知，距离误差、角误差与共面误差保持不变，即 $e_{12}^\dagger=e_{12}^\star,\ e_{14}^\dagger=e_{14}^\star,\ e_{\theta_i}^\dagger=e_{\theta_i}^\star,\ e_{cop}^\dagger=e_{cop}^\star$。若切平面基取为随 $Q$ 运输的 $E_d^\dagger=QE_d$，则法向误差分量也满足 $(e_n^\perp)^\dagger=(E_d^\dagger)^T(n_{142}^\dagger\times\bar n_d)=E_d^T(n_{142}^\star\times\bar n_d)=e_n^{\perp\star}$。

同理，速度误差满足 $\tilde v^\dagger=(I_4\otimes Q)\tilde v^\star$，而二阶积分闭环控制律中的每一项都只依赖于 $z_{ij},\ n_{234},\ n_{142},\ r_i$ 及各误差量，因此在刚体变换下满足等变性。于是有 $C^\dagger=C^\star(I_4\otimes Q^T)$，$(-B^\dagger)=(I_4\otimes Q)(-B^\star)$，从而 $G^\dagger=C^\dagger(-B^\dagger)=C^\star(-B^\star)=G^\star$，并且在最小坐标下有 $\hat G^\dagger=RG^\dagger S=\hat G^\star$，$\hat A^\dagger=T_Q\hat A^\star T_Q^{-1}$，其中 $T_Q=\operatorname{diag}(I_8,I_4\otimes Q)$。

**定理 2（任意同构目标方形的坐标无关局部指数稳定性）** 在定理 1 的同一二阶积分模型与增益条件下，对任意 $P^\dagger\in\mathcal S(\ell^\star,\bar n_d)$，若在 $P^\dagger$ 附近采用由 $Q$ 运输得到的局部最小内部误差坐标与切平面基，则相应最小内部误差系统的线性化矩阵 $\hat A^\dagger$ 与代表性目标方形处的 $\hat A^\star$ 共轭，从而具有相同谱。特别地，$P^\dagger$ 对应的最小内部误差平衡点也是局部指数稳定的。

**证明** 由上述刚体共轭关系可知 $\hat A^\dagger=T_Q\hat A^\star T_Q^{-1}$。由定理 1 知 $\hat A^\star$ 为 Hurwitz 矩阵，而 Hurwitz 性在相似变换下保持不变，故 $\hat A^\dagger$ 亦为 Hurwitz 矩阵。再对 $P^\dagger$ 附近的非线性系统重复定理 1 的论证，即得结论。证毕。

## 附录 A. 关键矩阵条目的复核路径

为便于审稿和答辩复核，本文关键矩阵的中间求导路径归纳如下。

1. 距离误差条目来自 $d_{ij}=\|p_j-p_i\|$ 的一阶微分，即 $D d_{ij}=z_{ij}^T\Pi_{ij}$。

2. 角误差条目来自 $\cos\theta_i=z_{ij}^Tz_{ik}$ 的链式求导，即 $D e_{\theta_i}=-(z_{ik}^TM_{ij}+z_{ij}^TM_{ik})$，并在代表性目标方形处使用 $\sin\theta_i^\star=1$。

3. 共面误差条目来自 $e_{cop}=z_{31}^Tn_{234}$ 的链式求导，即 $D e_{cop}=n_{234}^TM_{31}+z_{31}^TN_3(p)$。

4. 法向误差条目来自 $e_n^\perp=E_d^T(n_{142}\times\bar n_d)$ 的链式求导，即 $D e_n^\perp=-E_d^TS(\bar n_d)N_{124}(p)$。

5. $J_{sq}$ 与 $J_{raw}$ 由第 2 节中定义的映射 $\hat h(\xi)$ 与 $h_{raw}(\xi)$ 对局部坐标 $\xi$ 直接求导得到；$C^\star$ 则由本附录前四条在 $p^\star$ 处代入几何条件得到；$(-B^\star)$ 由二阶积分闭环控制律逐项在 $p^\star$ 处代入得到；$G^\star$ 与 $\hat G^\star$ 分别由矩阵乘积 $C^\star(-B^\star)$ 和 $RG^\star S$ 给出。

6. 配套的符号核验脚本位于 `analysis/derive_minimal_linearization.py`，可用于独立复核 $J_{sq}$、$J_{raw}$、$C^\star$、$(-B^\star)$、$G^\star$、$\hat G^\star$ 与特征多项式。
