## 3 Main Results

### 3.1 控制律设计 + 几何解释

控制律由五部分组成：

**角度项**
\[
u_1^{ang} = -k_\angle e_{\theta_1}(b_{12}+b_{14}), \quad
u_2^{ang} = -k_\angle e_{\theta_2}(b_{23}+b_{21}), \quad
u_3^{ang} = -k_\angle e_{\theta_3}(b_{34}+b_{32}), \quad
u_4^{ang} = -k_\angle e_{\theta_4}(b_{41}+b_{43})
\]
几何解释：每个角度项沿两个 bearing 的角平分线方向推拉，驱动内角恢复到 π/2。

**距离项**
\[
u_1^{dist} = k_d e_{12}b_{12} + k_d e_{14}b_{14}, \quad
u_2^{dist} = -k_d e_{12}b_{12}, \quad
u_3^{dist} = 0, \quad
u_4^{dist} = -k_d e_{14}b_{14}
\]
几何解释：沿 edge 方向伸缩，调节边长。只有 agent 1,2,4 参与距离控制。

**速度阻尼项**
\[
u_i^{vel} = -k_v(v_i - v_c^\star), \quad i \in \{1,2,3,4\}
\]

**共面项**
\[
u_3^{cop} = -k_{cop} e_{cop} n_{234}, \quad
u_1^{cop} = u_2^{cop} = u_4^{cop} = 0
\]
其中 \(e_{cop} = b_{31}^T n_{234}\)，\(n_{234} = (b_{32} \times b_{34}) / \|b_{32} \times b_{34}\|\)。
几何解释：沿 n_234 方向推 agent 3，将平面 {2,3,4} 倾斜以包含 agent 1。

**法向项**
\[
u_i^{nor} = k_n (e_n \times r_i), \quad i \in \{1,2,4\}, \quad u_3^{nor} = 0
\]
其中 \(e_n = n_{142} \times \bar n_d\)，\(r_i = p_i - (p_1+p_2+p_4)/3\)。
几何解释：\(e_n \times r_i\) 产生绕质心的力矩，旋转平面法向对齐目标，同时消解反射 ambiguity。

完整控制律：
\[
u_i = u_i^{ang} + u_i^{dist} + u_i^{vel} + u_i^{cop} + u_i^{nor}
\]

### 3.2 共面约束项的必要性 (Proposition)

> 仅依赖平面内角度约束不足以保证编队在三维机动中维持共面性。

**证明思路**：构造代表性目标方形，对 agent 3 施加离面扰动 h，计算角度和共面误差对 h 的敏感性：
- 角度误差：\(e_{\theta_3} = -h^2/(\ell^\star)^2 + O(h^4)\) → 二阶敏感
- 共面误差：\(e_{cop} = -h/(\sqrt{2}\ell^\star) + O(h^3)\) → 一阶敏感

结论：角度约束项对离面偏差缺乏足够的直接恢复力，必须显式引入共面恢复项。

### 3.3 名义稳定性证明 (Theorem 1)

**前提**：固定参考法向、固定目标分支、LESO 精确补偿（名义模型）。

**步骤**：
1. 构造最小内部误差坐标：证明 \(e_{\theta_3}\) 局部冗余（Lemma，见 Appendix）
2. 8 维最小误差 + 12 维速度误差 = 20 维闭环系统
3. 在代表性目标方形处线性化 → \(\hat A^\star\) 特征多项式分解
4. Hurwitz 判定：所有特征值严格位于 LHP
5. 局部指数稳定
6. 推广到同一固定分支上任意同构目标方形（刚体共轭）

### 3.4 LESO 观测器设计

对每个 agent i 的每个坐标轴 σ，设计三阶 LESO：
\[
\begin{aligned}
e_{i\sigma} &= z_{1,i\sigma} - p_{i\sigma} \\
\dot z_{1,i\sigma} &= z_{2,i\sigma} - \beta_{01} e_{i\sigma} \\
\dot z_{2,i\sigma} &= z_{3,i\sigma} + b_0 u_{i\sigma} - \beta_{02} e_{i\sigma} \\
\dot z_{3,i\sigma} &= -\beta_{03} e_{i\sigma}
\end{aligned}
\]
其中 \(\beta_{01}=3\omega_o,\ \beta_{02}=3\omega_o^2,\ \beta_{03}=\omega_o^3\)。

补偿输入：\(u_i = u_{0,i} - \hat D_i\)，其中 \(\hat D_i = z_{3,i}\) 估计总扰动（绳索张力 + 外部扰动 + 未建模动态）。

### 3.5 残差稳定性 (Proposition 1-3)

定义补偿残差 \(\varepsilon_i^{rope} = F_i^{rope}/m_i - \hat D_i\)。

- **Prop 1**：若 \(\varepsilon_i^{rope} \equiv 0\) → 真实外环退化为名义系统
- **Prop 2**：若初始分支正确、轨道在邻域内 → 迟滞逻辑不触发切换
- **Prop 3**：若 \(\varepsilon_i^{rope}\) 有界且小 → 局部 ISS / UUB；若 \(\varepsilon_i^{rope} \to 0\) → 局部渐近收敛
