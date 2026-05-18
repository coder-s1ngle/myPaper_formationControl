## 1 Introduction

### 1.1 研究背景

- 反无人机的社会需求越来越强烈，反无人机方式列举 -> 突出捕网方法的有效性(A concept for catching drones with a net carried by cooperative UAVs)
- 对于无人机编队的控制也已经做出了很多工作 -> 引出角度编队的优越性。
- 无人机编队可以有效提升反无人机效率。为了满足无人机编队反无人机的目的，需要无人机编队在三维空间中维持一个目标平面构型，因此不仅要求控制编队的法向量和尺度，还需要维持平面编队在三维空间中的共面特性。
- 特别的在捕网场景下，还要考虑柔性网对编队控制的影响


### 1.2 相关工作

- 目前关于角度编队的控制已经做了很多工作，但是平面编队控制往往需要给出智能体在同一平面下的假设，限制了编队在实际三维空间中的应用落地。
- 3D中的角度编队往往没有针对共面正方形的构型做出讨论(更多讨论多边体；henneberg的顶点添加方式，针对正方形四个agent共圆的情况不管用)
- 针对平面编队在3D空间中的控制也做了很多工作(leader-follower 相对位置 虚拟结构)，但是信息获取条件比较严格，当前只需要少量的距离信息和少量的通信，还有bearing测量，就可以在局部坐标系下完成编队控制，同时不需要局部坐标系对齐。


### 1.3 主要贡献


- 针对三维空间中的平面方形编队保持问题，指出仅依赖平面内角度约束不足以保证编队在机动过程中维持共面性，并据此引入显式共面恢复项，以抑制离面偏差并增强目标平面恢复能力。
- 面向捕网场景中网口平面朝向控制需求，提出一种无向平面法向控制方法。该方法仅依赖局部坐标系下的几何信息、少量距离测量和有限通信，无需全局坐标系对齐即可实现目标平面法向的局部对齐。
- 构造了四智能体平面方形编队的最小内部误差闭环系统,证明了名义闭环误差系统的局部指数稳定性.通过观测器补偿绳索张力以及外部扰动。当存在有界补偿误差的时候，证明闭环误差系统误差有界。


## 2 Problem formulation

### 2.1 Agent dynamics

- 给出四个 agent 在三维空间中的双积分模型：
  \[
  \dot p_i=v_i,\qquad \dot v_i=u_i,\qquad i\in\{1,2,3,4\}.
  \]
### 2.2 绳索简化模型
- 定义顶点集
  \[
  \mathcal V=\{1,2,3,4\}.
  \]
- 记 \(\mathcal N_i^{rope}\) 为与 agent \(i\) 存在绳索连接的邻居集合。
- 实际捕网场景中还存在弹性网张力。对与 agent \(i\) 相连的绳段 \(ij\) 定义
  \[
  \ell_{ij}=\|p_j-p_i\|,\qquad
  \dot \ell_{ij}=z_{ij}^\top(v_j-v_i),\qquad
  z_{ij}=\frac{p_j-p_i}{\|p_j-p_i\|},
  \]
  并采用弹簧-阻尼型张力模型
  \[
  T_{ij}=\big[k_r(\ell_{ij}-\ell_{ij}^0)+c_r\dot\ell_{ij}\big]_+
  \]
  其中 \([\,\cdot\,]_+=\max(\cdot,0)\)。
- 进一步给出作用在 agent \(i\) 上的绳索张力合力，通过环形绳索连接来描述简化的弹性网模型
  \[
  F_i^{rope}=\sum_{j\in\mathcal N_i^{rope}} T_{ij} z_{ij},
  \]

<!-- ### 2.2 Available measurements and desired planar square constraints

- 每个 agent 都能够在自身局部坐标系下获得控制律所需的 bearing 信息。
- 与尺度调节项相关的必要距离信息可由指定 agent 测得，并通过有限通信发送给相关邻居，以支持距离控制项及法向控制中相关几何量的构造。
- 期望公共速度和期望平面法向参考可用，或其在各 agent 局部坐标系下的表示可获得。
- 再引入用于控制设计的约束集合描述：
  - 顶点集 \(V=\{1,2,3,4\}\)
  - 角约束集合 \(A\)
  - 距离约束边集 \(E_d\)
- 定义目标网口构型对应的约束集合，例如：
  - 角约束集合 \(A=\{(2,1,4),(3,2,1),(4,3,2),(1,4,3)\}\)
  - 距离约束边集 \(E_d=\{(1,2),(1,4)\}\)
- 这些约束只刻画平面内的方形形状与尺度，还不足以完整描述三维空间中的目标平面网口。
- 定义目标平面方形集合 \(\mathcal S(\ell^\star,n_d)\)：除满足方形边长和角度约束外，还要求四点共面，并满足平面法向与目标法向在无向意义下一致。 -->

### 2.3 镜像反转问题说明并引出无向法向控制

- 在当前捕网场景下，更关心的是网口平面的无向朝向，而不是唯一有向刚体姿态。
- 因此目标法向按无向量处理，即 \(n_d\) 与 \(-n_d\) 代表同一目标平面朝向类。
- 说明本文不进一步区分 镜像 对应的两种法向符号，而是研究无向目标平面意义下的目标集合恢复。
- 引出无向法向分支选择：实际控制中从 \(\{n_d,-n_d\}\) 中选择当前局部一致的一支，记为 \(\bar n_d\)。

### 2.4 基于无向法向控制给出捕网场景下的期望构型描述

- 设 \(\ell^\star>0\) 为期望边长，\(n_d\in\mathbb R^3\) 为期望平面法向，且 \(\|n_d\|=1\)。
- 定义
  \[
  \mathcal S(\ell^\star,n_d)
  \]
  为目标平面方形构型集合,其满足：
  - 四个顶点构成方形，即四个内角均为 \(\pi/2\)；
  - 构型尺度由期望边长 \(\ell^\star\) 给定；
  - 四点共面；
  - 构型所在平面的法向与 \(n_d\) 一致，但按无向法向意义理解，即 \(n_d\) 与 \(-n_d\) 视为同一目标朝向类。
- 对于本文考虑的四智能体网口，平面内形状可由角约束集合
  \[
  A=\{(2,1,4),(3,2,1),(4,3,2),(1,4,3)\}
  \]
  描述，尺度可由距离约束边集
  \[
  E_d=\{(1,2),(1,4)\}
  \]
  给定。
- 还需进一步加入共面条件与无向法向条件。

### 2.5 Assumptions

#### Assumption 1
- 每个 agent 都能够在自身局部坐标系下获得控制律所需的 bearing 信息。

#### Assumption 2
- 与尺度调节项相关的必要距离信息可由指定 agent 测得，并通过有限通信发送给相关邻居，用以支持距离控制项及法向控制中相关几何量的构造。

#### Assumption 3
- 期望公共速度和期望平面法向参考可用，或其在各 agent 局部坐标系下的表示可获得。

#### Assumption 4
- 系统初始构型位于目标构型附近的充分小邻域内；期望距离大于 0，参考三角形非退化(不存在共线情况)。

### 2.6 几何量及误差表达方式

- 定义 bearing 向量、距离与相对位置：
  \[
  q_{ij}=p_j-p_i,\qquad d_{ij}=\|q_{ij}\|,\qquad z_{ij}=\frac{q_{ij}}{\|q_{ij}\|}.
  \]
- 定义无向内角
  \[
  \measuredangle jik=\arccos\left(z_{ij}^\top z_{ik}\right)\in[0,\pi].
  \]
- 共面项测量方式。
-  \(n_{142}\)测量方式。
- 定义无向法向误差 \(e_n\) 或其局部切平面分量 \(e_n^\perp\)。

### 2.7 Problem statement

#### Problem 1

- 设计分布式控制律，使得当前四智能体平面编队能够在3D空间中控制法向量还有距离缩放，并且保证其稳定性。
- 具体控制目标：
  - 角度误差收敛；
  - 距离/尺度误差收敛；
  - 共面误差收敛；
  - 无向法向误差收敛；
  - 速度误差收敛。


## 3 Main Results

### 3.1 控制律设计 + 几何解释

- 给出角度项、距离项、共面项、法向项、速度阻尼项设计。
- 明确各项的物理作用分工。
- 配合几何解释：
  - 角度项：bearings 沿角平分线方向推拉，恢复直角
  - 距离项：沿 edge 方向伸缩，调节边长
  - 速度阻尼项：跟踪公共速度
  - 共面项：沿 n_234 方向推 agent 3，将平面 {2,3,4} 倾斜以包含 agent 1
  - 法向项：叉乘 e_n × r_i 产生绕质心的力矩，旋转平面法向对齐目标

### 3.2 共面约束项的必要性 (Proposition) 

- 数学上说明角度约束项对离面误差不敏感，因此需要显式引入共面恢复项。
- 考虑代表性目标方形
  \[
  p_1^\star=[0,\ell^\star,0]^T,\quad
  p_2^\star=[0,0,0]^T,\quad
  p_3^\star=[\ell^\star,0,0]^T,\quad
  p_4^\star=[\ell^\star,\ell^\star,0]^T
  \]
  并构造单参数离面扰动
  \[
  p_3(h)=[\ell^\star,0,h]^T,\qquad h\in\mathbb R,\ |h|\ll 1,
  \]
  其余三个点保持不变。
- 先显式计算四个角。由
  \[
  p_2-p_1=[0,-\ell^\star,0]^T,\qquad p_4-p_1=[\ell^\star,0,0]^T
  \]
  可得
  \[
  \theta_1(h)=\frac{\pi}{2}.
  \]
  同理，由
  \[
  (p_3-p_2)^T(p_1-p_2)=0,\qquad
  (p_1-p_4)^T(p_3-p_4)=0
  \]
  可得
  \[
  \theta_2(h)=\frac{\pi}{2},\qquad
  \theta_4(h)=\frac{\pi}{2}.
  \]
- 对顶点 3，有
  \[
  p_4-p_3=[0,\ell^\star,-h]^T,\qquad
  p_2-p_3=[-\ell^\star,0,-h]^T,
  \]
  从而
  \[
  \cos\theta_3(h)=
  \frac{(p_4-p_3)^T(p_2-p_3)}
       {\|p_4-p_3\|\,\|p_2-p_3\|}
  =\frac{h^2}{(\ell^\star)^2+h^2}.
  \]
- 在 \(h=0\) 附近展开，得到
  \[
  \frac{h^2}{(\ell^\star)^2+h^2}
  =\frac{h^2}{(\ell^\star)^2}+O(h^4),
  \]
  以及
  \[
  \theta_3(h)=\arccos\!\left(\frac{h^2}{(\ell^\star)^2+h^2}\right)
  =\frac{\pi}{2}-\frac{h^2}{(\ell^\star)^2}+O(h^4).
  \]
  因而角误差满足
  \[
  e_{\theta_1}=e_{\theta_2}=e_{\theta_4}=0,\qquad
  e_{\theta_3}=-\frac{h^2}{(\ell^\star)^2}+O(h^4).
  \]
- 这一步要明确指出：对该代表性离面扰动，角度误差对离面变量 \(h\) 仅呈二阶敏感，而不是一阶敏感。
- 然后计算共面误差。若定义
  \[
  e_{cop}=z_{31}^T n_{234},
  \]
  其中
  \[
  z_{31}=\frac{p_1-p_3}{\|p_1-p_3\|}
  =\frac{[-\ell^\star,\ell^\star,-h]^T}{\sqrt{2(\ell^\star)^2+h^2}},
  \]
  而
  \[
  z_{32}=\frac{[-\ell^\star,0,-h]^T}{\sqrt{(\ell^\star)^2+h^2}},\qquad
  z_{34}=\frac{[0,\ell^\star,-h]^T}{\sqrt{(\ell^\star)^2+h^2}},
  \]
  则
  \[
  n_{234}
  =\frac{z_{32}\times z_{34}}{\|z_{32}\times z_{34}\|}
  =\frac{[h,-h,-\ell^\star]^T}{\sqrt{(\ell^\star)^2+2h^2}}.
  \]
- 于是
  \[
  e_{cop}(h)
  =\frac{[-\ell^\star,\ell^\star,-h]\cdot[h,-h,-\ell^\star]}
         {\sqrt{2(\ell^\star)^2+h^2}\sqrt{(\ell^\star)^2+2h^2}}
  =\frac{-\ell^\star h}{\sqrt{2(\ell^\star)^2+h^2}\sqrt{(\ell^\star)^2+2h^2}}.
  \]
  在 \(h=0\) 附近展开得到
  \[
  e_{cop}(h)
  =-\frac{1}{\sqrt2\,\ell^\star}h+O(h^3).
  \]
- 共面误差对离面位移 \(h\) 是一阶敏感，而角度误差仅二阶敏感。
- 单独依赖角度项时，离面模态缺乏足够强的直接恢复作用；
- 为保证平面编队在三维空间中的共面恢复与保持，需要显式引入共面约束项。

### 3.3 名义稳定性证明 (Theorem 1)

- 当前约束存在冗余：构造最小内部误差坐标（Lemma，放在 Appendix）
- 证明原始 9 维误差中 e_θ3 是冗余的（J_raw 秩=8）
- 构造 8 维最小内部误差 + 12 维速度误差 = 20 维闭环系统
- Hurwitz 判定 → 名义系统局部指数稳定
- 推广到固定分支上任意同构目标方形

### 3.4 LESO 观测器设计

- 给出 LESO 观测器设计，补偿实际捕网场景中的绳索张力与其他匹配扰动
- 三阶 LESO → 估计位置、速度、总扰动
- 补偿输入 u_i = u_{0,i} - \hat D_i
- 使实际外环在局部上逼近名义双积分外环

### 3.5 残差稳定性 (Proposition 1-3)

- 存在补偿残差 ε_i^{rope}
- 若 ε_i^{rope} ≡ 0 → 真实外环退化为名义系统（Prop 1）
- 若初始分支正确、轨道停留在邻域内 → 无切换（Prop 2）
- 若 ε_i^{rope} 小且有界 → 局部 ISS / 一致最终有界（Prop 3）
- 若 ε_i^{rope}(t) → 0 → 局部渐近收敛

### 4 Simulations
- 普通编队控制(去除绳索&LESO) 对比 去除共面约束项 -> 证明共面项的必要性。
- 加上绳索和LESO，给出捕捉case

### 5 Conclusion and discussion



def plot_3d_trajectory(agent_list, target_history=None, dt_sim=1/500,
                       t_snapshots=None):
    """
    Publication-quality 3-D trajectory plot for formation-control paper.

    Features
    --------
    - UAV trajectories (solid)
    - Formation-shape snapshots drawn as semi-transparent filled quads
      with edge wireframes
    - Time labels auto-positioned relative to formation centroid
    - Orthographic projection, clean pane-less style, balanced aspect

    Parameters
    ----------
    agent_list : list[Agent]
        Four QUAV agents with ``history_pos`` populated (NED frame).
    dt_sim : float
        Simulation time step (1 / 500 s by default).
    t_snapshots : list[float] or None
        Times at which to draw formation snapshots.  When None a
        reasonable default is chosen automatically.
    """
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    # ── colour palette (colourblind-friendly) ──
    c_uav = ["#2166AC", "#D6604D", "#4DAF4A", "#9467BD"]
    c_edge = "#D62728"
    c_face = "#CCCCCC"

    # ── NED → ENU ──
    pos_enu = [ned2enu_pos(np.array(a.history_pos)) for a in agent_list]

    # ── raw data ranges ──
    all_x = np.concatenate([p[:, 0] for p in pos_enu])
    all_y = np.concatenate([p[:, 1] for p in pos_enu])
    all_z = np.concatenate([p[:, 2] for p in pos_enu])

    x_min_data, x_max_data = np.min(all_x), np.max(all_x)
    y_min_data, y_max_data = np.min(all_y), np.max(all_y)
    z_min_data, z_max_data = np.min(all_z), np.max(all_z)

    rx = x_max_data - x_min_data
    ry = y_max_data - y_min_data
    rz = z_max_data - z_min_data

    # ── balance extreme aspect ratios ──
    r_max = max(rx, ry, rz, 1e-6)
    min_frac = 0.25
    rxv = max(rx, min_frac * r_max)
    ryv = max(ry, min_frac * r_max)
    rzv = 0.6 * max(rz, min_frac * r_max)
    mid_x = 0.5 * (x_min_data + x_max_data)
    mid_y = 0.5 * (y_min_data + y_max_data)
    mid_z = 0.5 * (z_min_data + z_max_data)

    xlim_lo = mid_x - 0.5 * rxv
    xlim_hi = mid_x + 0.5 * rxv
    ylim_lo = mid_y - 0.5 * ryv
    ylim_hi = mid_y + 0.5 * ryv
    zlim_lo = mid_z - 0.5 * rzv
    zlim_hi = mid_z + 0.5 * rzv

    # ── margins (percentage of balanced range) ──
    mx = 0.06 * rxv
    my = 0.06 * ryv
    mz = 0.10 * rzv

    # ── figure ──
    fig = plt.figure(figsize=(6.5, 4.0), dpi=300)
    ax = fig.add_subplot(111, projection='3d')
    ax.set_proj_type('ortho')
    ax.view_init(elev=15, azim=25)

    # ── box aspect from balanced ranges ──
    ax.set_box_aspect((rxv, ryv, rzv))

    # ── UAV trajectories ──
    for idx, (p, c) in enumerate(zip(pos_enu, c_uav)):
        ax.plot(p[:, 0], p[:, 1], p[:, 2],
                color=c, linewidth=0.5, linestyle='--', dashes=(4, 3),
                label=f'UAV {idx+1}')

    # ── default snapshot times ──
    if t_snapshots is None:
        T_end = len(pos_enu[0]) * dt_sim
        t_snapshots = np.linspace(0, T_end, 6)

    # ── formation snapshots ──
    for si, ts in enumerate(t_snapshots):
        idx_t = min(int(round(ts / dt_sim)), len(pos_enu[0]) - 1)
        verts = np.array([p[idx_t] for p in pos_enu])          # (4, 3)

        # filled quad
        poly = Poly3DCollection([list(verts)],
                                facecolor=c_face, edgecolor='none',
                                alpha=0.25, zorder=2)
        ax.add_collection3d(poly)

        # edge wireframe (perimeter)
        for (i, j) in [(0, 1), (1, 2), (2, 3), (3, 0)]:
            ax.plot([verts[i, 0], verts[j, 0]],
                    [verts[i, 1], verts[j, 1]],
                    [verts[i, 2], verts[j, 2]],
                    color=c_edge, linewidth=1.0, zorder=3)

        # net grid (井字形 thin lines representing lightweight capture net)
        v0, v1, v2, v3 = verts
        for frac in (1/3, 2/3):
            # lines parallel to edges 0-1 / 3-2 (one direction)
            p_a = v0 + frac * (v1 - v0)  # point on edge 0→1
            p_b = v3 + frac * (v2 - v3)  # point on edge 3→2
            ax.plot([p_a[0], p_b[0]], [p_a[1], p_b[1]], [p_a[2], p_b[2]],
                    color="#000000", linewidth=0.4, alpha=0.6, zorder=2)
            # lines parallel to edges 0-3 / 1-2 (other direction)
            p_c = v0 + frac * (v3 - v0)  # point on edge 0→3
            p_d = v1 + frac * (v2 - v1)  # point on edge 1→2
            ax.plot([p_c[0], p_d[0]], [p_c[1], p_d[1]], [p_c[2], p_d[2]],
                    color="#000000", linewidth=0.4, alpha=0.6, zorder=2)

        # time label — below formation, staggered left/right
        cx, cy, cz = verts[:, 0].mean(), verts[:, 1].mean(), verts[:, 2].mean()
        sign_x = -1.0 if si % 2 == 0 else +1.0
        label_dx = sign_x * 0.08 * rxv
        label_dy = 0.02 * ryv
        label_dz = 0.7 * rzv  #-0.15
        ax.text(cx + label_dx, cy + label_dy, cz + label_dz,
                rf'$t={ts:.0f}\,\mathrm{{s}}$',
                fontsize=5.0, color='#333333',
                ha='center', va='top', zorder=5)

    # ── limits ──
    ax.set_xlim(xlim_lo - mx, xlim_hi + mx)
    ax.set_ylim(ylim_lo - my, ylim_hi + my)
    ax.set_zlim(zlim_lo - mz, zlim_hi + mz)

    # ── axis labels ──
    ax.set_xlabel('X[m]', fontsize=6, labelpad=-3)
    ax.set_ylabel('Y[m]', fontsize=6, labelpad=2)
    ax.set_zlabel('Z[m]', fontsize=6, labelpad=-12)

    # ── ticks (per-axis lengths compensate 3-D projection foreshortening) ──
    ax.tick_params(labelsize=6, pad=2, direction='out', length=6)
    ax.xaxis.set_tick_params(length=1,pad = -2)
    ax.yaxis.set_tick_params(length=3,pad = -1)
    ax.zaxis.set_tick_params(length=100,pad = -4.5)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
    ax.zaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))

    # ── clean style ──
    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis._axinfo['grid'].update(color='#cccccc', alpha=0.4, linewidth=0.28)
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.set_visible(False)
    for spine in [ax.xaxis.line, ax.yaxis.line, ax.zaxis.line]:
        spine.set_linewidth(0.6)
        spine.set_color('black')

    # ── legend ──
    ax.legend(fontsize=4, ncol=4, loc='upper left',
              bbox_to_anchor=(0.01, 0.75), framealpha=0.85,
              edgecolor='gray', handlelength=1.0, handletextpad=0.3,
              columnspacing=0.6, borderpad=0.2, labelspacing=0.15)

    # ── layout: manual margins to keep axis labels visible ──
    fig.subplots_adjust(left=0.15, right=0.85, bottom=0.15, top=0.85)

    # ── save ──
    save_path = os.path.join(output_dir, 'trajectory_3d.pdf')
    fig.savefig(save_path, dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)
