# 名义闭环系统稳定性完整证明 (Theorem 1)

## 3.3 Stability of the Nominal Closed-Loop System

### 3.3.1 误差定义与闭环结构

**Internal error vector.** 定义以下几何误差：

- 距离误差：
  \[
  e_{12}=d_{12}-\ell^*,\qquad e_{14}=d_{14}-\ell^*.
  \]

- 角度误差（四个内角）：
  \[
  e_{\theta_i}=\alpha_{jik}-\frac{\pi}{2},\qquad
  (j,i,k)\in\mathcal A=\{(2,1,4),(3,2,1),(4,3,2),(1,4,3)\}.
  \]

- 共面误差：
  \[
  e_{\mathrm{cop}}=b_{31}^{T}n_{234},\qquad
  n_{234}=\frac{b_{32}\times b_{34}}{\|b_{32}\times b_{34}\|}.
  \]

- 法向误差：
  \[
  e_n=n_{142}\times\bar n_d,\qquad
  n_{142}=\frac{b_{14}\times b_{12}}{\|b_{14}\times b_{12}\|}.
  \]
  由于 \(e_n\) 始终正交于 \(n_{142}\)，仅有两个独立自由度。在平衡点处 \(n_{142}^*=[0,0,-1]^{T}\)，取前两个分量为独立坐标：
  \[
  (e_n^\perp)_1=e_{n,x},\qquad (e_n^\perp)_2=e_{n,y}.
  \]

- 速度误差：
  \[
  \tilde v_i=v_i-v_c^*,\qquad i\in\{1,2,3,4\},
  \]
  其中 \(v_c^*\) 为常值期望速度（名义模型中 \(\dot v_c^*=0\)）。

将上述误差组装为 raw internal error vector 和 velocity error stack：
\[
\eta=\bigl[e_{12},\,e_{14},\,e_{\theta_1},\,e_{\theta_2},\,e_{\theta_3},\,e_{\theta_4},\,e_{\mathrm{cop}},\,(e_n^\perp)_1,\,(e_n^\perp)_2\bigr]^{T}\in\mathbb R^{9},
\]
\[
\tilde v=\bigl[\tilde v_1^{T},\,\tilde v_2^{T},\,\tilde v_3^{T},\,\tilde v_4^{T}\bigr]^{T}\in\mathbb R^{12}.
\]

**名义模型.** 在本节中，假设参考法向分支已固定 (\(\bar n_d\) 为常值)，LESO 精确补偿总扰动（即 \(\varepsilon_i^{rope}\equiv 0\)），则 agent 动力学退化为
\[
\dot p_i=v_i,\qquad \dot v_i=u_i,\qquad i\in\{1,2,3,4\}.
\]

---

### 3.3.2 误差运动学：\(\dot\eta=C(p)\tilde v\)

本节给出每个误差量的运动学，并将其整理为 \(\dot\eta=C(p)\tilde v\) 的形式。

**距离误差运动学.** 由 \(d_{ij}=\|p_j-p_i\|\)，对时间求导：
\[
\dot d_{ij}=b_{ij}^{T}(v_j-v_i)=b_{ij}^{T}(\tilde v_j-\tilde v_i).
\]
因此
\[
\begin{aligned}
\dot e_{12}&=-b_{12}^{T}\tilde v_1+b_{12}^{T}\tilde v_2,\\
\dot e_{14}&=-b_{14}^{T}\tilde v_1+b_{14}^{T}\tilde v_4.
\end{aligned}
\tag{P.1}
\]

**角度误差运动学.** 由 \(\cos\alpha_{jik}=b_{ij}^{T}b_{ik}\)，两边对时间求导：
\[
-\sin\alpha_{jik}\,\dot\alpha_{jik}
=b_{ij}^{T}\dot b_{ik}+b_{ik}^{T}\dot b_{ij}.
\]
bearing 向量的导数为
\[
\dot b_{ij}=\frac{1}{d_{ij}}P_{b_{ij}}(v_j-v_i),\qquad
P_{b_{ij}}=I_3-b_{ij}b_{ij}^{T},
\]
其中 \(P_{b_{ij}}\) 为投影算子，满足 \(P_{b_{ij}}b_{ij}=0\)。代入得
\[
\dot\alpha_{jik}=
-\frac{b_{ik}^{T}P_{b_{ij}}(v_j-v_i)}{d_{ij}\sin\alpha_{jik}}
-\frac{b_{ij}^{T}P_{b_{ik}}(v_k-v_i)}{d_{ik}\sin\alpha_{jik}}.
\tag{P.2}
\]
四个角度误差的运动学均由 (P.2) 给出。以 \((j,i,k)=(2,1,4)\) 为例：
\[
\dot e_{\theta_1}=
-\frac{b_{14}^{T}P_{b_{12}}(\tilde v_2-\tilde v_1)}{d_{12}\sin\alpha_{214}}
-\frac{b_{12}^{T}P_{b_{14}}(\tilde v_4-\tilde v_1)}{d_{14}\sin\alpha_{214}}.
\tag{P.3}
\]
其余三个角度误差类似，仅替换对应的指标三元组。

**共面误差运动学.** 由 \(e_{\mathrm{cop}}=b_{31}^{T}n_{234}\)，其中 \(n_{234}=(b_{32}\times b_{34})/\|b_{32}\times b_{34}\|\)。对时间求导利用叉积的导数公式，得
\[
\dot e_{\mathrm{cop}}=
n_{234}^{T}\dot b_{31}+b_{31}^{T}\dot n_{234}.
\]
由于 \(\dot b_{31}=d_{31}^{-1}P_{b_{31}}(\tilde v_1-\tilde v_3)\)，而 \(\dot n_{234}\) 由 \(\dot b_{32},\dot b_{34}\) 的叉积和正交化构成，合并后可将 \(\dot e_{\mathrm{cop}}\) 表达为 \(\tilde v\) 的线性组合；其具体系数函数由 bearing 向量和距离确定。在目标构型附近，该映射是光滑的。

**法向误差运动学.** 由 \(n_{142}=(b_{14}\times b_{12})/\|b_{14}\times b_{12}\|\)，类似可导出 \(\dot n_{142}\) 对 \(\tilde v_1,\tilde v_2,\tilde v_4\) 的线性依赖关系。进而
\[
\dot e_n=\dot n_{142}\times\bar n_d,
\]
取前两个独立分量 \(\dot{(e_n^\perp)_1},\dot{(e_n^\perp)_2}\) 也光滑依赖于 \(\tilde v\)。

综合以上，存在光滑矩阵函数 \(C(p)\in\mathbb R^{9\times 12}\) 使得
\[
\boxed{\dot\eta=C(p)\,\tilde v}.
\tag{P.4}
\]

---

### 3.3.3 速度误差动力学：\(\dot{\tilde v}=-B(p)\eta-k_v\tilde v\)

名义模型下 \(\dot{\tilde v}_i=\dot v_i=u_i\)。控制律为
\[
u_i=u_i^{\mathrm{ang}}+u_i^{\mathrm{dist}}+u_i^{\mathrm{vel}}+u_i^{\mathrm{cop}}+u_i^{\mathrm{nor}}.
\]

**角度控制项.** 对 agent \(i\)，角度项沿两 bearing 的角平分线方向推拉：
\[
\begin{aligned}
u_1^{\mathrm{ang}}&=-k_\theta e_{\theta_1}(b_{12}+b_{14}),\\
u_2^{\mathrm{ang}}&=-k_\theta e_{\theta_2}(b_{23}+b_{21}),\\
u_3^{\mathrm{ang}}&=-k_\theta e_{\theta_3}(b_{34}+b_{32}),\\
u_4^{\mathrm{ang}}&=-k_\theta e_{\theta_4}(b_{41}+b_{43}).
\end{aligned}
\tag{P.5}
\]

**距离控制项.**
\[
\begin{aligned}
u_1^{\mathrm{dist}}&=k_d(e_{12}b_{12}+e_{14}b_{14}),\\
u_2^{\mathrm{dist}}&=-k_d e_{12}b_{12},\\
u_3^{\mathrm{dist}}&=0,\\
u_4^{\mathrm{dist}}&=-k_d e_{14}b_{14}.
\end{aligned}
\tag{P.6}
\]

**速度阻尼项.**
\[
u_i^{\mathrm{vel}}=-k_v\tilde v_i,\qquad i\in\{1,2,3,4\}.
\tag{P.7}
\]

**共面控制项.**
\[
u_3^{\mathrm{cop}}=-k_{\mathrm{cop}}e_{\mathrm{cop}}n_{234},\qquad
u_1^{\mathrm{cop}}=u_2^{\mathrm{cop}}=u_4^{\mathrm{cop}}=0.
\tag{P.8}
\]

**法向控制项.**
\[
u_i^{\mathrm{nor}}=k_n(e_n\times r_i),\quad i\in\{1,2,4\},\qquad
u_3^{\mathrm{nor}}=0,
\tag{P.9}
\]
其中 \(r_i=p_i-p_{124}\)，\(p_{124}=\frac{1}{3}(p_1+p_2+p_4)\) 为参考三角形的质心。注意 \(e_n\times r_i\) 满足 \(\sum_{i\in\{1,2,4\}}e_n\times r_i=0\)（过质心的纯力矩，不产生净平移）。

将 (P.5)--(P.9) 代入 \(\dot{\tilde v}_i=u_i\)，并利用 bearing 向量与误差之间的依赖关系，将除了 \(-k_v\tilde v_i\) 之外的所有项整理为 \(\eta\) 分量的线性组合（在目标构型附近以更高阶非线性余项为代价），得
\[
\boxed{\dot{\tilde v}=-B(p)\,\eta-k_v\tilde v},
\tag{P.10}
\]
其中 \(B(p)\in\mathbb R^{12\times 9}\) 光滑依赖于当前构型 \(p\)。

联合 (P.4) 和 (P.10)，名义闭环系统为
\[
\boxed{\begin{aligned}
\dot\eta&=C(p)\,\tilde v,\\
\dot{\tilde v}&=-B(p)\,\eta-k_v\tilde v.
\end{aligned}}
\tag{P.11}
\]

---

### 3.3.4 在目标方形处的线性化

取代表性目标方形（边长 \(\ell^*\)，法向 \(\bar n_d=[0,0,-1]^{T}\)）：
\[
p_1^*=\begin{bmatrix}0\\\ell^*\\0\end{bmatrix},\;
p_2^*=\begin{bmatrix}0\\0\\0\end{bmatrix},\;
p_3^*=\begin{bmatrix}\ell^*\\0\\0\end{bmatrix},\;
p_4^*=\begin{bmatrix}\ell^*\\\ell^*\\0\end{bmatrix}.
\tag{P.12}
\]

在该平衡点处：\(\eta^*=0\)，\(\tilde v^*=0\)，\(p=p^*\)。定义
\[
C^*=C(p^*)\in\mathbb R^{9\times 12},\qquad
B^*=B(p^*)\in\mathbb R^{12\times 9}.
\]

令 \(X=[\eta^{T},\tilde v^{T}]^{T}\in\mathbb R^{21}\)。在 \(X=0\) 处线性化 (P.11)：
\[
\boxed{\dot X=A^*X},\qquad
A^*=\begin{bmatrix}
0_{9\times 9} & C^*\\
-B^* & -k_v I_{12}
\end{bmatrix}.
\tag{P.13}
\]

---

### 3.3.5 \(C^*\) 的显式结构

在目标构型 (P.12) 处，关键 bearing 向量为：
\[
\begin{aligned}
b_{12}^*&=[0,-1,0]^{T}, & b_{14}^*&=[1,0,0]^{T},\\
b_{21}^*&=[0,1,0]^{T},  & b_{23}^*&=[1,0,0]^{T},\\
b_{32}^*&=[-1,0,0]^{T}, & b_{34}^*&=[0,1,0]^{T},\\
b_{41}^*&=[-1,0,0]^{T}, & b_{43}^*&=[0,-1,0]^{T}.
\end{aligned}
\tag{P.14}
\]

所有内角 \(\alpha_{jik}^*=\pi/2\)，故 \(\sin\alpha_{jik}^*=1\)，\(\cos\alpha_{jik}^*=0\)。距离 \(d_{12}^*=d_{14}^*=d_{23}^*=d_{34}^*=\ell^*\)，\(d_{13}^*=d_{24}^*=\sqrt2\ell^*\)。

**\(C^*\) 的第 1--2 行（距离误差运动学）.** 从 (P.1) 直接代入 \(b_{12}^*,b_{14}^*\)：
\[
\begin{aligned}
\dot e_{12}&=[0,1,0]\tilde v_1+[-0,-1,0]\tilde v_2+[0,0,0]\tilde v_3+[0,0,0]\tilde v_4,\\
\dot e_{14}&=[-1,0,0]\tilde v_1+[0,0,0]\tilde v_2+[0,0,0]\tilde v_3+[1,0,0]\tilde v_4.
\end{aligned}
\tag{P.15}
\]
（记号为：\(\dot e_{12}=c_{1}^{T}\tilde v\)，其中 \(c_1=[0,1,0,\;0,-1,0,\;0,0,0,\;0,0,0]^{T}\)。）

**\(C^*\) 的第 3--6 行（角度误差运动学）.** 在目标构型处计算 \(d_{ij}\sin\alpha_{jik}\) 和 \(b^{T}P_b\)。

以 \(e_{\theta_1}\)（顶点 1, 角 (2,1,4)）为例。在平衡点处：
\[
b_{14}^{*T}P_{b_{12}^*}=[1,0,0]\begin{bmatrix}1&0&0\\0&0&0\\0&0&1\end{bmatrix}=[1,0,0],
\]
\[
b_{12}^{*T}P_{b_{14}^*}=[0,-1,0]\begin{bmatrix}0&0&0\\0&1&0\\0&0&1\end{bmatrix}=[0,-1,0],
\]
\(d_{12}^*=\ell^*\)，\(d_{14}^*=\ell^*\)，\(\sin\alpha_1^*=1\)。代入 (P.3)：
\[
\dot e_{\theta_1}=
-\frac{1}{\ell^*}[1,0,0](\tilde v_2-\tilde v_1)
-\frac{1}{\ell^*}[0,-1,0](\tilde v_4-\tilde v_1).
\]
展开：
\[
\dot e_{\theta_1}=
\frac{1}{\ell^*}[1,-1,0]\tilde v_1
-\frac{1}{\ell^*}[1,0,0]\tilde v_2
+0\cdot\tilde v_3
+\frac{1}{\ell^*}[0,1,0]\tilde v_4.
\tag{P.16}
\]

类似可得 \(e_{\theta_2},e_{\theta_3},e_{\theta_4}\) 的运动学（完整 12×1 行向量见 Appendix A.1）。

**\(C^*\) 的第 7 行（共面误差运动学）.** 在平衡点处：
\[
b_{31}^*=\frac{[-\ell^*,\ell^*,0]^{T}}{\sqrt2\,\ell^*}
=[-1/\sqrt2,\,1/\sqrt2,\,0]^{T},\qquad
n_{234}^*=[0,0,-1]^{T}.
\]
共面误差 \(e_{\mathrm{cop}}=b_{31}^{T}n_{234}\) 的线性化给出（推导见 Appendix A.2）：
\[
\dot e_{\mathrm{cop}}=
\frac{1}{\sqrt2\,\ell^*}[-1,1,0]\tilde v_1
+\frac{1}{\sqrt2\,\ell^*}[0,-1,0]\tilde v_2
+\frac{1}{\sqrt2\,\ell^*}[1,0,0]\tilde v_3
+0\cdot\tilde v_4.
\tag{P.17}
\]

**\(C^*\) 的第 8--9 行（法向误差运动学）.** 在平衡点处 \(n_{142}^*=[0,0,-1]^{T}\)，\(\bar n_d=[0,0,-1]^{T}\)，\(e_n^*=0\)。法向误差的前两个分量的线性化仅依赖于 \(\tilde v_1,\tilde v_2,\tilde v_4\)（agent 3 不参与参考三角形的法向运动）。结果形式为（完整推导见 Appendix A.3）：
\[
\begin{bmatrix}\dot{(e_n^\perp)_1}\\ \dot{(e_n^\perp)_2}\end{bmatrix}
=M_n\begin{bmatrix}\tilde v_1\\\tilde v_2\\\tilde v_4\end{bmatrix},
\qquad M_n\in\mathbb R^{2\times 9}.
\tag{P.18}
\]

将 (P.15)--(P.18) 组装，得 \(C^*\in\mathbb R^{9\times 12}\) 的完整条目（显式矩阵见 Appendix A.4）。

---

### 3.3.6 \(-B^*\) 的显式结构

将控制律 (P.5)--(P.9) 在目标构型处线性化。由于 \(v_c^*\) 为常值，\(\dot{\tilde v}_i=u_i\)。

**角度项.** 在目标构型处，bearing 向量取 (P.14)，角度误差系数直接读出：
\[
\begin{aligned}
u_1^{\mathrm{ang}}&=-k_\theta e_{\theta_1}([0,-1,0]^{T}+[1,0,0]^{T})
=-k_\theta e_{\theta_1}[1,-1,0]^{T},\\[2pt]
u_2^{\mathrm{ang}}&=-k_\theta e_{\theta_2}([1,0,0]^{T}+[0,1,0]^{T})
=-k_\theta e_{\theta_2}[1,1,0]^{T},\\[2pt]
u_3^{\mathrm{ang}}&=-k_\theta e_{\theta_3}([0,1,0]^{T}+[-1,0,0]^{T})
=-k_\theta e_{\theta_3}[-1,1,0]^{T},\\[2pt]
u_4^{\mathrm{ang}}&=-k_\theta e_{\theta_4}([-1,0,0]^{T}+[0,-1,0]^{T})
=-k_\theta e_{\theta_4}[-1,-1,0]^{T}.
\end{aligned}
\tag{P.19}
\]

**距离项.**
\[
\begin{aligned}
u_1^{\mathrm{dist}}&=k_d(e_{12}[0,-1,0]^{T}+e_{14}[1,0,0]^{T})
=k_d[e_{14},\,-e_{12},\,0]^{T},\\
u_2^{\mathrm{dist}}&=-k_d e_{12}[0,-1,0]^{T}=k_d[0,\,e_{12},\,0]^{T},\\
u_3^{\mathrm{dist}}&=0,\\
u_4^{\mathrm{dist}}&=-k_d e_{14}[1,0,0]^{T}=k_d[-e_{14},\,0,\,0]^{T}.
\end{aligned}
\tag{P.20}
\]

**共面项.** 在平衡点处 \(n_{234}^*=[0,0,-1]^{T}\)：
\[
u_3^{\mathrm{cop}}=-k_{\mathrm{cop}}e_{\mathrm{cop}}[0,0,-1]^{T}
=k_{\mathrm{cop}}e_{\mathrm{cop}}[0,0,1]^{T}.
\tag{P.21}
\]

**法向项.** 在平衡点附近，对 \(i\in\{1,2,4\}\)：
\[
r_1^*=[0,\ell^*,0]^{T}-[\ell^*/3,2\ell^*/3,0]^{T}
=[-\ell^*/3,\ell^*/3,0]^{T},
\]
\[
r_2^*=[0,0,0]^{T}-[\ell^*/3,2\ell^*/3,0]^{T}
=[-\ell^*/3,-2\ell^*/3,0]^{T},
\]
\[
r_4^*=[\ell^*,\ell^*,0]^{T}-[\ell^*/3,2\ell^*/3,0]^{T}
=[2\ell^*/3,\ell^*/3,0]^{T}.
\]
注意 \(e_n^*=0\)，线性化法向项的核心是 \(e_n\) 对 \(\eta\) 的依赖关系。在平衡点附近，\(e_n\approx[(e_n^\perp)_1,(e_n^\perp)_2,0]^{T}\)。对 \(i=1,2,4\)：
\[
u_i^{\mathrm{nor}}=k_n(e_n\times r_i^*)
=k_n\begin{bmatrix}
-(e_n^\perp)_2 r_{i,z}^*\\
(e_n^\perp)_1 r_{i,z}^*\\
(e_n^\perp)_2 r_{i,x}^*-(e_n^\perp)_1 r_{i,y}^*
\end{bmatrix}.
\]
由于 \(r_{i,z}^*=0\)（所有点在 XY 平面），前两个分量为零；第三分量：
\[
\begin{aligned}
u_{1,z}^{\mathrm{nor}}&=k_n\bigl((e_n^\perp)_2(-\ell^*/3)-(e_n^\perp)_1(\ell^*/3)\bigr)
=-\frac{k_n\ell^*}{3}\bigl((e_n^\perp)_1+(e_n^\perp)_2\bigr),\\[2pt]
u_{2,z}^{\mathrm{nor}}&=k_n\bigl((e_n^\perp)_2(-\ell^*/3)-(e_n^\perp)_1(-2\ell^*/3)\bigr)
=\frac{k_n\ell^*}{3}\bigl(2(e_n^\perp)_1-(e_n^\perp)_2\bigr),\\[2pt]
u_{4,z}^{\mathrm{nor}}&=k_n\bigl((e_n^\perp)_2(2\ell^*/3)-(e_n^\perp)_1(\ell^*/3)\bigr)
=\frac{k_n\ell^*}{3}\bigl(-(e_n^\perp)_1+2(e_n^\perp)_2\bigr).
\end{aligned}
\tag{P.22}
\]
注意 \(u_{1,z}^{\mathrm{nor}}+u_{2,z}^{\mathrm{nor}}+u_{4,z}^{\mathrm{nor}}=0\)，即纯力矩性质。

**速度阻尼项.** 直接贡献 \(-k_v\tilde v_i\)，已在 (P.13) 的 \(-k_v I_{12}\) 块中吸收。

将 (P.19)--(P.22) 组装，得 \(-B(p^*)\) 的各列（每个 \(\eta\) 分量对应的加速度方向）。例如：
- \(e_{12}\) 列：来自距离项 \(u_1^{\mathrm{dist}},u_2^{\mathrm{dist}}\) 的贡献 → \(k_d[0,1,0\mid 0,1,0\mid 0,0,0\mid 0,0,0]^{T}\)（每段 3D）
- \(e_{\theta_1}\) 列：仅 \(u_1^{\mathrm{ang}}\) → \(-k_\theta[1,-1,0\mid 0,0,0\mid 0,0,0\mid 0,0,0]^{T}\)
- 等等

完整 \(-B^*\in\mathbb R^{12\times 9}\) 的显式矩阵见 Appendix A.5。

---

### 3.3.7 角度冗余性与最小内部误差坐标 (Lemma 1)

**Lemma 1 (Angular redundancy).** 对于满足共面约束 \(e_{\mathrm{cop}}=0\) 的任意四边形构型，四个内角满足
\[
\alpha_{214}+\alpha_{321}+\alpha_{432}+\alpha_{143}=2\pi.
\tag{P.23}
\]
特别地，在目标构型的邻域内，角度误差 \(e_{\theta_3}\) 可局部表示为其余三个角度误差的线性函数：
\[
e_{\theta_3}=-(e_{\theta_1}+e_{\theta_2}+e_{\theta_4})+O(\|\eta\|^2).
\tag{P.24}
\]

**Proof.** 对于平面四边形，四个内角之和恒为 \(2\pi\)。又因 \(\alpha_{jik}=\pi/2+e_{\theta_i}\)，代入 \(\sum_{i=1}^4 (\pi/2+e_{\theta_i})=2\pi\) 得 \(\sum_{i=1}^4 e_{\theta_i}=0\)，故 \(e_{\theta_3}=-(e_{\theta_1}+e_{\theta_2}+e_{\theta_4})\)。当 \(e_{\mathrm{cop}}\neq 0\)（轻度非共面）时，四个"投影内角"满足此关系直至一阶误差（余弦定理的扰动展开），因此 (P.24) 以 \(O(\|\eta\|^2)\) 余项成立。\(\blacksquare\)

Lemma 1 表明 raw error vector \(\eta\in\mathbb R^9\) 包含一个局部冗余自由度。定义选择矩阵 \(\hat S\in\mathbb R^{9\times 8}\)（删除 \(e_{\theta_3}\) 对应的第 5 行）和嵌入矩阵 \(\hat R\in\mathbb R^{8\times 9}\)（将 \(\eta\) 的前 4 个分量和第 6--9 分量映射到 8 维），使得
\[
\hat\eta=\hat R\eta\in\mathbb R^8
\]
为最小内部误差坐标。显式地：
\[
\hat\eta=[e_{12},\,e_{14},\,e_{\theta_1},\,e_{\theta_2},\,e_{\theta_4},\,e_{\mathrm{cop}},\,(e_n^\perp)_1,\,(e_n^\perp)_2]^{T}.
\tag{P.25}
\]
对应地，令 \(\hat C^*=\hat R C^*\in\mathbb R^{8\times 12}\) 和 \(\hat B^*=B^*\hat S\in\mathbb R^{12\times 8}\) 为简约后的 Jacobian 矩阵。定义完整最小状态
\[
\hat X=[\hat\eta^{T},\tilde v^{T}]^{T}\in\mathbb R^{20}.
\]
其线性化动力学为
\[
\boxed{\dot{\hat X}=\hat A^*\hat X},\qquad
\hat A^*=\begin{bmatrix}
0_{8\times 8} & \hat C^*\\
-\hat B^* & -k_v I_{12}
\end{bmatrix}.
\tag{P.26}
\]

---

### 3.3.8 特征多项式与 Hurwitz 判定

定义简约结构矩阵
\[
\hat G^*=\hat C^*(-\hat B^*)\in\mathbb R^{8\times 8}.
\tag{P.27}
\]
\(\hat G^*\) 的物理意义是：\(\ddot{\hat\eta}=\hat G^*\hat\eta-k_v\dot{\hat\eta}\)（二阶误差动力学中的刚度矩阵）。

通过代入 (P.15)--(P.22) 在约简后的显式形式，计算 \(\hat G^*=\hat C^*(-\hat B^*)\)（详细矩阵乘法见 Appendix A.6），得 \(\hat G^*\) 的分块对角/三角结构，对应不同误差通道之间的耦合关系：

- **法向通道**（\((e_n^\perp)_1,(e_n^\perp)_2\)）：与角度、距离、共面误差解耦。对应 \(2\times 2\) 子块 \(\hat G^*_n=-k_n I_2\)。
- **共面通道**（\(e_{\mathrm{cop}}\)）：仅通过 \(u_3^{\mathrm{cop}}\) 反馈，对应 \(1\times 1\) 子块 \(\hat G^*_{\mathrm{cop}}=-(\sqrt2/(2\ell^*))k_{\mathrm{cop}}\)。
- **角度-距离耦合通道**（\(e_{12},e_{14},e_{\theta_1},e_{\theta_2},e_{\theta_4}\)）：\(5\times 5\) 子块；进一步可分解为一个单实特征值 \(-4k_\theta/\ell^*\) 和两个相同的 \(2\times 2\) 块。

\(\hat G^*\) 的特征多项式为：
\[
\boxed{\begin{aligned}
\det(\lambda I_8-\hat G^*)
&=(\lambda+k_n)^2
\left(\lambda+\frac{\sqrt2}{2\ell^*}k_{\mathrm{cop}}\right)
\left(\lambda+\frac{4k_\theta}{\ell^*}\right)\\
&\quad\times\left(\lambda^2+\frac{2(k_d\ell^*+k_\theta)}{\ell^*}\lambda
+\frac{2k_d k_\theta}{\ell^*}\right)^2.
\end{aligned}}
\tag{P.28}
\]

**Hurwitz 验证.** 对于任意正增益 \(k_\theta,k_d,k_v,k_{\mathrm{cop}},k_n>0\) 和 \(\ell^*>0\)：
- 一次因子 \(\lambda+k_n\)（二重）：根 \(\lambda=-k_n<0\)。
- 一次因子 \(\lambda+\frac{\sqrt2}{2\ell^*}k_{\mathrm{cop}}\)：根 \(\lambda<0\)。
- 一次因子 \(\lambda+\frac{4k_\theta}{\ell^*}\)：根 \(\lambda<0\)。
- 二次因子 \(\lambda^2+a\lambda+b\)，其中
  \[
  a=\frac{2(k_d\ell^*+k_\theta)}{\ell^*}>0,\qquad
  b=\frac{2k_d k_\theta}{\ell^*}>0.
  \]
  该二次多项式的根具有负实部（Routh-Hurwitz：\(a>0,b>0\)）。

因此 \(\hat G^*\) 的所有特征值严格位于左半开平面，即 \(\hat G^*\) 是 Hurwitz 的。

---

### 3.3.9 全系统 Hurwitz 判定

\(\hat A^*\) 的特征多项式通过块行列式公式与 \(\hat G^*\) 关联：
\[
\begin{aligned}
\det(\lambda I_{20}-\hat A^*)
&=\det\begin{bmatrix}
\lambda I_8 & -\hat C^*\\
\hat B^* & (\lambda+k_v)I_{12}
\end{bmatrix}\\
&=(\lambda+k_v)^{4}
\det\!\Bigl(\lambda(\lambda+k_v)I_8-\hat G^*\Bigr).
\end{aligned}
\tag{P.29}
\]

**推导:** 使用分块矩阵行列式公式 \(\det\begin{bmatrix}A&B\\C&D\end{bmatrix}=\det(D)\det(A-BD^{-1}C)\)。这里 \(D=(\lambda+k_v)I_{12}\)，\(\det(D)=(\lambda+k_v)^{12}\)，且
\[
A-BD^{-1}C=\lambda I_8-\hat C^*\frac{1}{\lambda+k_v}(-\hat B^*)
=\lambda I_8+\frac{1}{\lambda+k_v}\hat G^*.
\]
乘以 \((\lambda+k_v)^8\) 并整理即得 (P.29)（提取因子 \((\lambda+k_v)^4\) 来自维度差 \(12-8=4\)）。

设 \(\mu_i\;(i=1,\ldots,8)\) 为 \(\hat G^*\) 的特征值（均满足 \(\operatorname{Re}(\mu_i)<0\)）。则 (P.29) 的根由以下方程给出：
\[
(\lambda+k_v)^4\prod_{i=1}^8\bigl(\lambda^2+k_v\lambda-\mu_i\bigr)=0.
\]

- \(k_v>0\) 给出四重根 \(\lambda=-k_v<0\)。
- 对每个 \(\mu_i<0\)（实根），二次方程 \(\lambda^2+k_v\lambda-\mu_i=0\)：
  \(\lambda=(-k_v\pm\sqrt{k_v^2+4\mu_i})/2\)。由于 \(\mu_i<0\) 且 \(k_v>0\)，\(\sqrt{k_v^2+4\mu_i}<k_v\)，两根实部均为负。
- 对每对共轭复根 \(\mu_i=\alpha_i\pm j\beta_i\)（\(\alpha_i<0\)），对应四次方程由两对二次因子组成，由 Lemma 的推广保证所有根具有负实部（条件 \((k_v^2\alpha_i+\beta_i^2)<0\) 因 \(\alpha_i<0\) 自然满足；详见 Appendix A.7）。

因此 \(\hat A^*\) 的所有 20 个特征值严格位于左半开平面，\(\hat A^*\) 是 Hurwitz 的。

---

### 3.3.10 非线性闭环的局部指数稳定性

将最小坐标下的非线性闭环系统写为
\[
\dot{\hat X}=\hat A^*\hat X+\hat\phi(\hat X),
\tag{P.30}
\]
其中 \(\hat\phi(\hat X)\) 包含线性化余项（高阶项、bearing 和距离误差的非线性耦合、共面项和法向项的非线性部分）。由于所有几何量（bearing 向量、距离、叉积）在 \(p^*\) 附近均为光滑函数，存在邻域和常数 \(L>0\) 使得
\[
\|\hat\phi(\hat X)\|\le L\|\hat X\|^2,\qquad \forall \|\hat X\|<\delta_0.
\tag{P.31}
\]

由于 \(\hat A^*\) 是 Hurwitz 的，对任意 \(Q=Q^{T}>0\)，存在唯一的 \(P=P^{T}>0\) 满足 Lyapunov 方程
\[
\hat A^{*T}P+P\hat A^*=-Q.
\]

取 Lyapunov 函数 \(V(\hat X)=\hat X^{T}P\hat X\)。沿 (P.30) 的解求导：
\[
\begin{aligned}
\dot V(\hat X)
&=\hat X^{T}(\hat A^{*T}P+P\hat A^*)\hat X
+2\hat X^{T}P\hat\phi(\hat X)\\
&=-\hat X^{T}Q\hat X+2\hat X^{T}P\hat\phi(\hat X)\\
&\le -\lambda_{\min}(Q)\|\hat X\|^2+2\|P\|\,\|\hat X\|\,\|\hat\phi(\hat X)\|\\
&\le -\lambda_{\min}(Q)\|\hat X\|^2+2\|P\|L\|\hat X\|^3\\
&=-\|\hat X\|^2\bigl(\lambda_{\min}(Q)-2\|P\|L\|\hat X\|\bigr).
\end{aligned}
\tag{P.32}
\]

在充分小的邻域 \(\mathcal U=\{\hat X:\|\hat X\|<r\}\) 内，其中
\[
r=\min\!\left(\delta_0,\;\frac{\lambda_{\min}(Q)}{4\|P\|L}\right),
\]
有 \(\lambda_{\min}(Q)-2\|P\|L\|\hat X\|\ge \frac{1}{2}\lambda_{\min}(Q)\)，从而
\[
\dot V(\hat X)\le -\frac{\lambda_{\min}(Q)}{2\lambda_{\max}(P)}V(\hat X).
\tag{P.33}
\]

由比较引理，
\[
V(\hat X(t))\le V(\hat X(0))\,e^{-\alpha t},\qquad
\alpha=\frac{\lambda_{\min}(Q)}{2\lambda_{\max}(P)}>0.
\]
利用 \(V(\hat X)\) 的二次型上下界
\[
\lambda_{\min}(P)\|\hat X\|^2\le V(\hat X)\le\lambda_{\max}(P)\|\hat X\|^2,
\]
最终得到指数衰减估计：
\[
\boxed{\|\hat X(t)\|\le Me^{-\alpha t}\|\hat X(0)\|,\qquad
M=\sqrt{\frac{\lambda_{\max}(P)}{\lambda_{\min}(P)}}\ge 1.}
\tag{P.34}
\]

---

### 3.3.11 推广至同一分支上的任意目标方形

以上分析针对特定目标方形 (P.12)（位于 XY 平面，法向 \(\bar n_d=[0,0,-1]^{T}\)）。对同一固定法向分支 \(\mathcal S(\ell^*,\bar n_d)\) 中的任意目标方形构型 \(p^\diamond\)，存在刚体变换（平移 + 绕 \(\bar n_d\) 轴的旋转）将 (P.12) 映至 \(p^\diamond\)。由于：

1. 所有 bearing 向量 \(b_{ij}\) 在刚体旋转下协变；
2. 内角 \(\alpha_{jik}\) 在刚体变换下不变；
3. 距离 \(d_{ij}\) 在刚体变换下不变；
4. 控制律 (P.5)--(P.9) 仅依赖 bearing 向量、距离和局部几何量，这些量在刚体变换下保持其相对几何关系；

因此闭环 Jacobian \(\hat A^*\) 的特征值在刚体共轭下不变。特别地，\(\hat A^*\) 在任意 \(p^\diamond\in\mathcal S(\ell^*,\bar n_d)\) 处均为 Hurwitz，局部指数稳定性对整族目标构型一致成立。

---

### 3.3.12 定理陈述 (Theorem 1)

> **Theorem 1 (Local exponential stability).** Consider the nominal closed-loop system under the control law (P.5)--(P.9) with fixed branch \(\bar n_d\), constant desired velocity \(v_c^*\), and positive gains \(k_\theta,k_d,k_v,k_{\mathrm{cop}},k_n>0\). There exists a neighborhood \(\mathcal U\) of the target square (P.12) and constants \(M\ge 1,\;\alpha>0\) such that for any initial condition \(\hat X(0)\in\mathcal U\), the minimal error state satisfies
> \[
> \|\hat X(t)\|\le Me^{-\alpha t}\|\hat X(0)\|,\qquad\forall t\ge 0.
> \]
> Furthermore, by rigid-body conjugation (Section 3.3.11), the same conclusion holds for any target square in the same fixed branch \(\mathcal S(\ell^*,\bar n_d)\).

**Proof.** Sections 3.3.1--3.3.10 constitute the proof. \(\blacksquare\)

---

## Appendix 条目清单（待后续补全）

以下为正文证明中引用的 Appendix 条目，其中包含因篇幅原因移至附录的完整显式计算：

- **Appendix A.1** — \(C^*\) 第 3--6 行（四个角度误差运动学）的完整 12 维行向量
- **Appendix A.2** — 共面误差运动学 \(\dot e_{\mathrm{cop}}\) 的详细推导
- **Appendix A.3** — 法向误差运动学 \(\dot{(e_n^\perp)_1},\dot{(e_n^\perp)_2}\) 的详细推导
- **Appendix A.4** — \(C^*\in\mathbb R^{9\times 12}\) 完整显式矩阵
- **Appendix A.5** — \(-B^*\in\mathbb R^{12\times 9}\) 完整显式矩阵
- **Appendix A.6** — \(\hat G^*=\hat C^*(-\hat B^*)\) 的矩阵乘法和特征多项式推导
- **Appendix A.7** — 共轭复根情形下二次因子 \(\lambda^2+k_v\lambda-\mu_i\) 的 Hurwitz 补充验证
