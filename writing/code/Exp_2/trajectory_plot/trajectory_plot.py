import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from matplotlib.ticker import MaxNLocator
from matplotlib.path import Path
import matplotlib
import os

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = ['Times New Roman']
matplotlib.rcParams['mathtext.fontset'] = 'stix'

output_dir = None


def set_output_dir(path):
    global output_dir
    output_dir = os.fspath(path)
    os.makedirs(output_dir, exist_ok=True)


def _crop_white_border_and_save(fig, save_path, margin_h=0.15, margin_v=0.15, dpi=300):
    """
    Save figure with controlled white borders as true vector PDF.

    1. Render to PNG with bbox_inches='tight' + safe pad → pixel-detect content
    2. Map pixel bounds back to figure coordinates (inches)
    3. Apply ``margin_h`` / ``margin_v`` ratio padding
    4. Save as vector PDF with computed bbox
    """
    import io
    from PIL import Image
    from matplotlib.transforms import Bbox

    safe_pad = 0.30

    # Step 1: render to PNG for pixel-level content detection
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, facecolor='white', edgecolor='white',
                bbox_inches='tight', pad_inches=safe_pad)
    buf.seek(0)

    img = Image.open(buf).convert('RGB')
    arr = np.array(img)
    h_px, w_px = arr.shape[:2]

    # Step 2: detect content bounding box
    mask = np.any(arr < 250, axis=2)
    rows = np.any(mask, axis=1)
    cols = np.any(mask, axis=0)

    if not rows.any() or not cols.any():
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight', pad_inches=0.10)
        buf.close()
        return

    rmin, rmax = np.where(rows)[0][[0, -1]]
    cmin, cmax = np.where(cols)[0][[0, -1]]

    # Step 3: map pixel bounds → figure coordinates (inches)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    tight_bbox = fig.get_tightbbox(renderer)

    px_to_in = 1.0 / dpi
    # PNG row 0 = top; figure y increases upward
    cx0 = tight_bbox.x0 - safe_pad + cmin * px_to_in
    cy0 = tight_bbox.y0 - safe_pad + (h_px - rmax) * px_to_in
    cx1 = tight_bbox.x0 - safe_pad + cmax * px_to_in
    cy1 = tight_bbox.y0 - safe_pad + (h_px - rmin) * px_to_in

    # Step 4: apply controlled margins
    cw = cx1 - cx0
    ch = cy1 - cy0

    final_bbox = Bbox.from_extents(
        cx0 - cw * margin_h,
        cy0 - ch * margin_v,
        cx1 + cw * margin_h,
        cy1 + ch * margin_v,
    )

    # Step 5: save true vector PDF
    fig.savefig(save_path, dpi=dpi, bbox_inches=final_bbox, pad_inches=0)
    buf.close()


# ==========================================
#         Coordinate Conversion Helpers
# ==========================================
def ned2enu_pos(pos_ned):
    """
    Convert position from NED (North-East-Down) to ENU (East-North-Up).
    pos_ned: (N, 3) or (3,) numpy array
    Returns: pos_enu with X=East, Y=North, Z=Up
    """
    if pos_ned.ndim == 1:
        # Single point: [x_ned, y_ned, z_ned] -> [y_ned, x_ned, -z_ned]
        return np.array([pos_ned[1], pos_ned[0], -pos_ned[2]])
    else:
        # Array of points
        pos_enu = np.zeros_like(pos_ned)
        pos_enu[:, 0] = pos_ned[:, 1]  # X_enu = Y_ned (East)
        pos_enu[:, 1] = pos_ned[:, 0]  # Y_enu = X_ned (North)
        pos_enu[:, 2] = -pos_ned[:, 2] # Z_enu = -Z_ned (Up)
        return pos_enu

def ned2enu_yaw(yaw_ned):
    """
    Convert Yaw from NED (0=North, Clockwise) to ENU (0=East, Counter-Clockwise).
    Formula: yaw_enu = pi/2 - yaw_ned
    """
    return np.pi/2 - yaw_ned

# ==========================================
#          Post-process Helpers
# ==========================================

def _crop_3d_figure(filepath, margin_px=4):
    """
    Auto-crop whitespace from a saved 3-D matplotlib figure.

    Tightly trims content, then pastes it onto a fresh white canvas
    with the requested margins.  This guarantees exact margin control
    even when ``bbox_inches='tight'`` already pushed content to the edge.

    Parameters
    ----------
    margin_px : int | tuple[int, int, int, int]
        Single int for uniform margin, or (top, bottom, left, right).
    """
    from PIL import Image
    import numpy as np

    if isinstance(margin_px, int):
        mt = mb = ml = mr = margin_px
    else:
        mt, mb, ml, mr = margin_px

    img = Image.open(filepath)
    arr = np.array(img)
    h, w = arr.shape[:2]

    if arr.ndim == 3 and arr.shape[2] == 4:
        rgb = arr[:, :, :3]
    else:
        rgb = arr

    white = np.all(rgb > 250, axis=2)

    # find tight content bounds
    top = 0
    while top < h and np.all(white[top, :]):
        top += 1

    bottom = h - 1
    while bottom > top and np.all(white[bottom, :]):
        bottom -= 1

    left = 0
    while left < w and np.all(white[:, left]):
        left += 1

    right = w - 1
    while right > left and np.all(white[:, right]):
        right -= 1

    # crop to content only
    content = img.crop((left, top, right + 1, bottom + 1))

    # paste onto fresh canvas with exact margins
    cw, ch = content.size
    new_w = cw + ml + mr
    new_h = ch + mt + mb
    canvas = Image.new('RGB', (new_w, new_h), (255, 255, 255))
    canvas.paste(content, (ml, mt))
    canvas.save(filepath)


# ==========================================
#               Plotting Functions
# ==========================================

def draw_uav_rectangle(ax, pos, yaw, uav_id="", size=0.16, color="blue", z_fixed=None):
    """
    绘制无人机矩形 + 编号
    Input pos and yaw should be in the plotting frame (ENU).
    """
    x, y, z = pos
    if z_fixed is not None:
        z = z_fixed

    # 本体矩形形状 (Body Frame)
    pts_body = np.array([
        [ size,  size],
        [ size, -size],
        [-size, -size],
        [-size,  size],
        [ size,  size]
    ])

    # Rotation matrix for ENU (Standard mathematical rotation)
    # R = [[cos, -sin], [sin, cos]]
    R = np.array([
        [ np.cos(yaw), -np.sin(yaw)],
        [ np.sin(yaw),  np.cos(yaw)]
    ])

    pts_world = (R @ pts_body.T).T
    pts_world[:,0] += x
    pts_world[:,1] += y

    # 绘制矩形
    if isinstance(ax, Axes3D):
        ax.plot(pts_world[:,0], pts_world[:,1], [z]*pts_world.shape[0], color=color, linewidth=2)
        # ax.text(x, y, z+0.1, uav_id, color="red", fontsize=12)
    else:
        ax.plot(pts_world[:,0], pts_world[:,1], color=color, linewidth=2)
        # ax.text(x, y, uav_id, color="red", fontsize=12)


def plot_xy_trajectory(agent_list):
    """
    绘制 XY 真实轨迹 + 预期轨迹 + 最终无人机矩形 (ENU Frame: East-North)
    同时保存为 xy_trajectory.png
    """
    plt.figure()
    ax = plt.gca()

    all_x, all_y = [], []   # 用来统计全局范围

    for agent in agent_list:
        if len(agent.history_pos) == 0:
            continue

        # Get NED data
        pos_ned = np.array(agent.history_pos)
        att_hist = np.array(agent.history_att) # [Roll, Pitch, Yaw] in NED
        yaw_ned = att_hist[:, 2]

        # Convert to ENU
        pos_enu = ned2enu_pos(pos_ned)
        yaw_enu = ned2enu_yaw(yaw_ned)

        # 收集范围
        all_x.extend(pos_enu[:, 0])
        all_y.extend(pos_enu[:, 1])

        # === 真实轨迹 (East, North) ===
        ax.plot(pos_enu[:, 0], pos_enu[:, 1], label=f"{agent.agent_id}_real")

        # === 预期轨迹 (Convert NED cmd to ENU) ===
        if hasattr(agent, "history_cmd") and len(agent.history_cmd) > 0:
            cmd_ned = np.array(agent.history_cmd)
            cmd_enu = ned2enu_pos(cmd_ned)
            ax.plot(cmd_enu[:, 0], cmd_enu[:, 1], '--', label=f"{agent.agent_id}_cmd")

        # 最终位置矩形
        draw_uav_rectangle(
            ax,
            pos_enu[-1],
            yaw_enu[-1],
            uav_id=agent.agent_id,
            size=0.16,
            color="black"
        )

    ax.set_xlabel("East [m]")
    ax.set_ylabel("North [m]")
    ax.set_title("2D Trajectory (ENU Frame)")
    ax.grid(True)
    ax.legend()

    # ======= 自适应坐标轴比例 =======
    if len(all_x) > 0 and len(all_y) > 0:
        xmin, xmax = min(all_x), max(all_x)
        ymin, ymax = min(all_y), max(all_y)

        x_range = max(xmax - xmin, 1e-6)
        y_range = max(ymax - ymin, 1e-6)

        # 给一点边距
        margin_x = 0.05 * x_range
        margin_y = 0.05 * y_range
        ax.set_xlim(xmin - margin_x, xmax + margin_x)
        ax.set_ylim(ymin - margin_y, ymax + margin_y)

        ratio = x_range / y_range

        # 如果两轴范围差不多，就保持等比例；否则用自适应
        if 1/3.0 <= ratio <= 3.0:
            ax.set_aspect("equal", adjustable="box")
        else:
            ax.set_aspect("auto")

    # === 保存图片 ===
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,"xy_trajectory.pdf"), dpi=300, bbox_inches="tight")
def _drone_path():
    """Return a simple quadcopter-shaped matplotlib Path marker."""
    r = 1.0       # arm length
    rb = 0.2      # body radius
    rr = 0.18     # rotor radius
    verts = []
    codes = []
    # four arms (lines from center to rotor positions) + rotors at corners
    rotors = [(r, r), (-r, r), (-r, -r), (r, -r)]
    for rx, ry in rotors:
        verts.extend([(rx, ry), (0, 0)])
        codes.extend([Path.MOVETO, Path.LINETO])
        # small circle for each rotor
        n = 10
        for i in range(n + 1):
            a = 2 * np.pi * i / n
            verts.append((rx + rr * np.cos(a), ry + rr * np.sin(a)))
            codes.append(Path.LINETO if i > 0 else Path.MOVETO)
    # center body
    n = 12
    for i in range(n + 1):
        a = 2 * np.pi * i / n
        verts.append((rb * np.cos(a), rb * np.sin(a)))
        codes.append(Path.LINETO if i > 0 else Path.MOVETO)
    return Path(verts, codes)


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
                label=f'QUAV {idx+1}')

    # ── default snapshot times ──
    if t_snapshots is None:
        T_end = len(pos_enu[0]) * dt_sim
        t_snapshots = np.linspace(0, T_end, 6)

    # ── formation snapshots ──
    for si, ts in enumerate(t_snapshots):
        idx_t = min(int(round(ts / dt_sim)), len(pos_enu[0]) - 1)
        verts = np.array([p[idx_t] for p in pos_enu])          # (4, 3)

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
                    color='#999999', linewidth=0.4, alpha=0.6, zorder=2)
            # lines parallel to edges 0-3 / 1-2 (other direction)
            p_c = v0 + frac * (v3 - v0)  # point on edge 0→3
            p_d = v1 + frac * (v2 - v1)  # point on edge 1→2
            ax.plot([p_c[0], p_d[0]], [p_c[1], p_d[1]], [p_c[2], p_d[2]],
                    color='#999999', linewidth=0.4, alpha=0.6, zorder=2)

        # captured target within net (t >= 65 s)
        if ts >= 65 and target_history is not None:
            net_ctr = np.mean(verts, axis=0)
            l_mid = 0.5 * (verts[0] + verts[1])
            r_mid = 0.5 * (verts[2] + verts[3])
            lo_mid = 0.5 * (verts[1] + verts[2])
            u_mid = 0.5 * (verts[0] + verts[3])
            ax_u = r_mid - l_mid; ax_v = u_mid - lo_mid
            hu = 0.5 * np.linalg.norm(ax_u); hv = 0.5 * np.linalg.norm(ax_v)
            eu = ax_u / (2 * hu + 1e-12); ev = ax_v / (2 * hv + 1e-12)
            cp = net_ctr + 0.25 * hu * eu - 0.15 * hv * ev
            sz = 2 if ts >= 85 else 4
            ax.scatter(*cp, s=sz, color='black', marker='.',
                       zorder=8)

        # time label — below formation, staggered left/right
        cx, cy, cz = verts[:, 0].mean(), verts[:, 1].mean(), verts[:, 2].mean()
        sign_x = -1.0 if si % 2 == 0 else +1.0
        label_dx = sign_x * 0.08 * rxv
        label_dy = 0.02 * ryv
        label_dz = 0.7 * rzv
        # per-snapshot nudges to avoid overlap
        if ts == 70:
            label_dx -= 0.35 * rxv
            label_dz -= 0.20 * rzv
        ax.text(cx + label_dx, cy + label_dy, cz + label_dz,
                rf'$t={ts:.0f}\,\mathrm{{s}}$',
                fontsize=5.0, color='#333333',
                ha='center', va='top', zorder=5)

    # ── target trajectory (drawn last to stay on top) ──
    if target_history is not None and len(target_history) > 0:
        tgt_enu = ned2enu_pos(np.array(target_history))
        ax.plot(tgt_enu[:, 0], tgt_enu[:, 1], tgt_enu[:, 2],
                color='black', linewidth=0.6, linestyle='--', dashes=(3, 2),
                label='Target', zorder=10)
        drone_marker = _drone_path()
        n_tgt = len(tgt_enu)
        # markers at start, middle, end — size decreasing start → end
        ax.scatter(*tgt_enu[0], s=18, color='black',
                   marker=drone_marker, zorder=10, facecolors='none')
        ax.scatter(*tgt_enu[len(tgt_enu) // 2], s=11, color='black',
                   marker=drone_marker, zorder=10, facecolors='none')
        ax.scatter(*tgt_enu[-1], s=7, color='black',
                   marker=drone_marker, zorder=20, facecolors='none')

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
        axis._axinfo['grid'].update(color='#cccccc', alpha=0.4, linewidth=0.05)
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

    # ── save with controlled white borders ──
    save_path = os.path.join(output_dir, 'trajectory_3d.pdf')
    _crop_white_border_and_save(fig, save_path, margin_h=0.05, margin_v=0.05, dpi=300)
    plt.close(fig)



def plot_paper_3d_trajectory(agent_list, target_history, t_sample=(0, 60, 80, 120)):
    """
    Publication-quality 3-D trajectory plot for the target-capture scenario.

    Parameters
    ----------
    agent_list : list[Agent]
        Four QUAV agents with history_pos, history_att, history_time populated.
    target_history : np.ndarray of shape (N, 3)
        Target positions in NED frame, sampled at 500 Hz.
    t_sample : tuple of float
        Times (seconds) at which to draw formation-shape snapshots.
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from matplotlib.colors import to_rgba

    # ── colour palette (colourblind-friendly, muted) ──
    c_uav = ["#2166AC", "#D6604D", "#4DAF4A", "#9467BD"]  # blue, red-orange, green, purple
    c_target = "#222222"

    # ── NED → ENU ──
    pos_enu = []
    for agent in agent_list:
        pos_enu.append(ned2enu_pos(np.array(agent.history_pos)))
    tgt_enu = ned2enu_pos(np.array(target_history))

    # ── figure setup ──
    fig = plt.figure(figsize=(5.2, 4.0), dpi=300)
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(elev=22, azim=-48)

    # ── UAV trajectories + ground shadows ──
    for idx, (p, c) in enumerate(zip(pos_enu, c_uav)):
        ax.plot(p[:, 0], p[:, 1], p[:, 2],
                color=c, linewidth=0.8, label=f'QUAV {idx+1}')
        # ground projection
        ax.plot(p[:, 0], p[:, 1], 0,
                color=to_rgba(c, 0.15), linewidth=0.4, zorder=0)

    # ── target trajectory + shadow ──
    ax.plot(tgt_enu[:, 0], tgt_enu[:, 1], tgt_enu[:, 2],
            color=c_target, linewidth=0.8, linestyle='--', dashes=(4, 3),
            label='Target')
    ax.plot(tgt_enu[:, 0], tgt_enu[:, 1], 0,
            color=to_rgba(c_target, 0.12), linewidth=0.4, zorder=0)

    # ── formation snapshots at t_sample ──
    dt = 1 / 500
    for ts in t_sample:
        idx = int(ts / dt)
        if idx >= len(pos_enu[0]):
            idx = len(pos_enu[0]) - 1
        verts = np.array([p[idx] for p in pos_enu])
        # close the quad
        verts_closed = np.vstack([verts, verts[0]])
        poly3d = [[tuple(v) for v in verts]]
        alpha_val = 0.35 if ts == 60 else 0.18
        face = Poly3DCollection(poly3d, alpha=alpha_val,
                                facecolor='lightgray', edgecolor='gray',
                                linewidth=0.5, zorder=2)
        ax.add_collection3d(face)

        # label the snapshot
        cx, cy, cz = verts[:, 0].mean(), verts[:, 1].mean(), verts[:, 2].mean()
        offset = np.array([0.3, -0.6, 0.4])
        ax.text(cx + offset[0], cy + offset[1], cz + offset[2],
                f'$t={ts}\\,\\mathrm{{s}}$',
                fontsize=7, color='gray', ha='left', va='bottom')

    # ── target position markers ──
    for ts in t_sample:
        idx = int(ts / dt)
        if idx >= len(tgt_enu):
            idx = len(tgt_enu) - 1
        tx, ty, tz = tgt_enu[idx]
        ax.scatter(tx, ty, tz, marker='*', color=c_target, s=40, zorder=5)

    # ── axes ──
    ax.set_xlabel('East / m', fontsize=9, labelpad=6)
    ax.set_ylabel('North / m', fontsize=9, labelpad=6)
    ax.set_zlabel('Height / m', fontsize=9, labelpad=6)
    ax.tick_params(labelsize=7, pad=2)
    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis._axinfo['grid'].update(color='#cccccc', alpha=0.4, linewidth=0.05)

    # ── legend ──
    ax.legend(fontsize=7, loc='upper left', ncol=2,
              framealpha=0.85, edgecolor='gray', handlelength=1.8)

    # ── set equal aspect on XY, tighten layout ──
    all_e = np.concatenate([p[:, 0] for p in pos_enu] + [tgt_enu[:, 0]])
    all_n = np.concatenate([p[:, 1] for p in pos_enu] + [tgt_enu[:, 1]])
    e_mid = 0.5 * (all_e.min() + all_e.max())
    n_mid = 0.5 * (all_n.min() + all_n.max())
    half = 0.55 * max(np.ptp(all_e), np.ptp(all_n))
    ax.set_xlim(e_mid - half, e_mid + half)
    ax.set_ylim(n_mid - half, n_mid + half)

    # z limits — keep tight around actual height range
    all_h = np.concatenate([p[:, 2] for p in pos_enu])
    h_range = np.ptp(all_h)
    h_mid = 0.5 * (all_h.min() + all_h.max())
    pad_z = max(0.3, 0.2 * h_range)
    ax.set_zlim(h_mid - pad_z, h_mid + pad_z)

    save_path = os.path.join(output_dir, 'trajectory_3d.pdf')
    _crop_white_border_and_save(fig, save_path, margin_h=0.10, margin_v=0.10, dpi=300)
    plt.close(fig)


def plot_attitude(agent_list):
    """
    以 subplot 方式绘制多个 UAV 的 roll/pitch/yaw 角度曲线 (NED Frame)
    3 行 1 列： roll, pitch, yaw
    同时保存为 attitude_angles.png
    """
    fig, axs = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

    titles = ["Roll (rad)", "Pitch (rad)", "Yaw (rad)"]

    for agent in agent_list:
        att_hist = np.array(agent.history_att)
        t = np.array(agent.history_time)

        if len(att_hist) == 0:
            continue

        # roll
        axs[0].plot(t, att_hist[:,0], label=agent.agent_id)
        # pitch
        axs[1].plot(t, att_hist[:,1], label=agent.agent_id)
        # yaw
        axs[2].plot(t, att_hist[:,2], label=agent.agent_id)

    for i in range(3):
        axs[i].set_ylabel(titles[i])
        axs[i].grid(True)
        axs[i].legend()

    axs[2].set_xlabel("Time [s]")

    fig.suptitle("UAV Attitude Angles (NED Frame)", fontsize=14)
    plt.tight_layout()

    # === 保存图片 ===
    fig.savefig(os.path.join(output_dir,"attitude_angles.pdf"), dpi=300, bbox_inches="tight")


def plot_z_over_time(agent_list):
    """
    绘制高度曲线 (Height = -Z_ned)
    同时保存为 z_height.png
    """
    plt.figure()
    ax = plt.gca()

    for agent in agent_list:
        if len(agent.history_pos) == 0:
            continue

        t = np.array(agent.history_time)
        pos_ned = np.array(agent.history_pos)
        
        # Height is -Z in NED
        height = -pos_ned[:, 2]

        ax.plot(t, height, label=f"{agent.agent_id}_height")

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Height [m]")
    ax.set_title("Height Over Time (ENU Z)")
    ax.grid(True)
    ax.legend()

    # === 保存图片 ===
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,"z_height.pdf"), dpi=300, bbox_inches="tight")


def make_gif(agent_list, filename="trajectory.gif", interval=20, sample_step=20):
    """
    2D 轨迹动图 (ENU Frame)
    """

    # ======== PRECOMPUTE GLOBAL AXIS RANGE (ENU) ========
    all_x = []
    all_y = []

    for agent in agent_list:
        pos_ned = np.array(agent.history_pos)
        if len(pos_ned) > 0:
            pos_enu = ned2enu_pos(pos_ned)
            all_x.extend(pos_enu[:, 0]) # East
            all_y.extend(pos_enu[:, 1]) # North

    if not all_x:
        print("No trajectory data.")
        return

    margin = 1.0
    xmin, xmax = min(all_x) - margin, max(all_x) + margin
    ymin, ymax = min(all_y) - margin, max(all_y) + margin

    size = 0.16
    # Body frame for visualization (will be rotated)
    rect_body = np.array([
        [ size,  size],
        [ size, -size],
        [-size, -size],
        [-size,  size],
        [ size,  size]
    ])

    max_len = max(len(a.history_pos) for a in agent_list)
    sampled = list(range(0, max_len, sample_step))
    if sampled[-1] != max_len - 1:
        sampled.append(max_len - 1)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_title("UAV Cluster Animation (ENU Frame)")
    ax.set_xlabel("East [m]")
    ax.set_ylabel("North [m]")
    ax.grid(True)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    # ======= 自适应/等比例选择 =======
    x_range = max(xmax - xmin, 1e-6)
    y_range = max(ymax - ymin, 1e-6)
    ratio = x_range / y_range
    if 1/3.0 <= ratio <= 3.0:
        ax.set_aspect("equal", adjustable="box")
    else:
        ax.set_aspect("auto")

    traj_lines = []
    uav_lines = []
    labels = []

    for agent in agent_list:
        line, = ax.plot([], [], linewidth=1.5, label=agent.agent_id)
        traj_lines.append(line)

        rect_plot, = ax.plot([], [], color="black", linewidth=2)
        uav_lines.append(rect_plot)

        txt = ax.text(0, 0, agent.agent_id, color="red", fontsize=12)
        labels.append(txt)

    ax.legend()

    def init():
        for line in traj_lines:
            line.set_data([], [])
        for rect in uav_lines:
            rect.set_data([], [])
        return traj_lines + uav_lines + labels

    def update(frame_idx):
        frame = sampled[frame_idx]

        for k, agent in enumerate(agent_list):
            pos_hist = agent.history_pos
            att_hist = agent.history_att

            if frame >= len(pos_hist):
                continue
            
            # Convert current frame to ENU
            cur_pos_ned = pos_hist[frame]
            cur_yaw_ned = att_hist[frame][2]
            
            cur_pos_enu = ned2enu_pos(cur_pos_ned)
            cur_yaw_enu = ned2enu_yaw(cur_yaw_ned)

            x, y = cur_pos_enu[0], cur_pos_enu[1]
            yaw = cur_yaw_enu

            # Convert trajectory history to ENU
            # Slicing up to frame
            hist_ned_segment = np.array(pos_hist[:frame+1])
            if len(hist_ned_segment) > 0:
                hist_enu_segment = ned2enu_pos(hist_ned_segment)
                traj_lines[k].set_data(hist_enu_segment[:, 0], hist_enu_segment[:, 1])
            else:
                traj_lines[k].set_data([], [])

            # Draw Drone Rect
            R = np.array([
                [np.cos(yaw), -np.sin(yaw)],
                [np.sin(yaw),  np.cos(yaw)]
            ])
            pts = (R @ rect_body.T).T
            pts[:, 0] += x
            pts[:, 1] += y
            uav_lines[k].set_data(pts[:, 0], pts[:, 1])

            labels[k].set_position((x, y))

        return traj_lines + uav_lines + labels

    anim = animation.FuncAnimation(
        fig, update, init_func=init,
        frames=len(sampled), interval=interval, blit=True
    )

    anim.save(os.path.join(output_dir, filename), writer=animation.PillowWriter(fps=1000 // interval))

    plt.close(fig)


def make_3d_gif(agent_list, filename="trajectory_3d.gif", interval=20, sample_step=20):
    """
    生成 3D 轨迹动图 (ENU Frame)
    """
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # ======== 1. 计算全局坐标范围 (ENU) ========
    all_x, all_y, all_z = [], [], []
    for agent in agent_list:
        pos_ned = np.array(agent.history_pos)
        if len(pos_ned) > 0:
            pos_enu = ned2enu_pos(pos_ned)
            all_x.extend(pos_enu[:, 0])
            all_y.extend(pos_enu[:, 1])
            all_z.extend(pos_enu[:, 2])

    if not all_x:
        print("No trajectory data found.")
        return

    margin = 0.5
    xmin, xmax = min(all_x) - margin, max(all_x) + margin
    ymin, ymax = min(all_y) - margin, max(all_y) + margin
    zmin, zmax = min(all_z) - margin, max(all_z) + margin
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_zlim(zmin, zmax)
    
    ax.set_xlabel("East [m]")
    ax.set_ylabel("North [m]")
    ax.set_zlabel("Height [m]")
    ax.set_title("UAV 3D Cluster Animation (ENU Frame)")

    # ======== 2. 初始化绘图元素 ========
    traj_lines = []
    uav_lines = []
    labels = []
    
    size = 0.16
    rect_body = np.array([
        [ size,  size],
        [ size, -size],
        [-size, -size],
        [-size,  size],
        [ size,  size]
    ])

    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

    for i, agent in enumerate(agent_list):
        c = colors[i % len(colors)]
        
        # 轨迹线
        line, = ax.plot([], [], [], linewidth=1.5, label=agent.agent_id, color=c)
        traj_lines.append(line)
        
        # 无人机本体 (用线框表示)
        uav_line, = ax.plot([], [], [], linewidth=2, color=c)
        uav_lines.append(uav_line)
        
        # 标签
        txt = ax.text(0, 0, 0, agent.agent_id, color='black')
        labels.append(txt)

    ax.legend()

    # 采样帧数
    max_len = max(len(a.history_pos) for a in agent_list)
    sampled = list(range(0, max_len, sample_step))
    if sampled[-1] != max_len - 1:
        sampled.append(max_len - 1)

    # ======== 3. 动画更新函数 ========
    def update(frame_idx):
        frame = sampled[frame_idx]
        
        for k, agent in enumerate(agent_list):
            pos_hist = agent.history_pos
            att_hist = agent.history_att
            
            if frame >= len(pos_hist):
                continue
                
            # --- Update Trajectory (ENU) ---
            hist_ned = np.array(pos_hist[:frame+1])
            if len(hist_ned) > 0:
                hist_enu = ned2enu_pos(hist_ned)
                xs = hist_enu[:, 0]
                ys = hist_enu[:, 1]
                zs = hist_enu[:, 2]
                
                traj_lines[k].set_data(xs, ys)
                traj_lines[k].set_3d_properties(zs)
            
            # --- Update Drone (ENU) ---
            cur_pos_ned = pos_hist[frame]
            cur_yaw_ned = att_hist[frame][2]
            
            cur_pos_enu = ned2enu_pos(cur_pos_ned)
            cur_yaw_enu = ned2enu_yaw(cur_yaw_ned)
            
            x, y, z = cur_pos_enu[0], cur_pos_enu[1], cur_pos_enu[2]
            yaw = cur_yaw_enu
            
            # Rotation matrix (ENU)
            R = np.array([
                [np.cos(yaw), -np.sin(yaw)],
                [np.sin(yaw),  np.cos(yaw)]
            ])
            
            # Calculate rect points in ENU
            pts = (R @ rect_body.T).T
            pts[:, 0] += x
            pts[:, 1] += y
            
            uav_lines[k].set_data(pts[:, 0], pts[:, 1])
            uav_lines[k].set_3d_properties([z] * 5) # Fixed height at current z
            
            # --- Update Label ---
            labels[k].set_position((x, y))
            labels[k].set_3d_properties(z + 0.2)

        return traj_lines + uav_lines + labels

    # ======== 4. 生成并保存 ========
    print("Generating 3D GIF, please wait...")
    anim = animation.FuncAnimation(
        fig, update, frames=len(sampled), interval=interval, blit=False
    )
    
    save_path = os.path.join(output_dir, filename)
    anim.save(save_path, writer=animation.PillowWriter(fps=1000 // interval))
    print(f"3D GIF saved to {save_path}")
    plt.close(fig)


# ==========================================
#          Data Analysis Helpers
# ==========================================

def make_time_axis(history_list, dt):
    """
    history_list: dict[str, list]
    dt: float
    """
    N = len(next(iter(history_list.values())))
    return np.arange(N) * dt


def plot_formation_edge_err(history_list, dt):
    t = make_time_axis(history_list, dt)

    plt.figure()
    for name, data in history_list.items():
        plt.plot(t, data, label=name)

    plt.xlabel("Time [s]")
    plt.ylabel("Edge distance [m]")
    plt.title("Formation Edge Distance Evolution")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,"formation_edge_history.pdf"), dpi=300, bbox_inches="tight")

def plot_formation_angle_err(history_list, dt):
    t = make_time_axis(history_list, dt)

    plt.figure()
    for name, data in history_list.items():
        plt.plot(t, data, label=name)

    plt.xlabel("Time [s]")
    plt.ylabel("Angle error [rad]")
    plt.title("Formation Angle Error Evolution")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,"formation_angle_err_history.pdf"), dpi=300, bbox_inches="tight")

def plot_formation_tension(history_list, dt):
    t = make_time_axis(history_list, dt)

    plt.figure()
    for name, data in history_list.items():
        plt.plot(t, data, label=name)

    plt.xlabel("Time [s]")
    plt.ylabel("Tension [N]")
    plt.title("Rope Tension Evolution")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,"formation_tension_history.pdf"), dpi=300, bbox_inches="tight")



def plot_Fd(agent_list):
    plt.figure()
    for agent in agent_list:
        t = np.array(agent.history_time)
        plt.plot(t, agent.history_Fd, label=agent.agent_id)
    plt.xlabel("Time [s]")
    plt.ylabel("Fd on agent [N]")
    plt.title("Fd on agent Evolution")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,"Fd_history.pdf"), dpi=300, bbox_inches="tight")


def plot_u(agent_list):
    fig,axs = plt.subplots(3,1,figsize = (10,10),sharex = True)

    titles = ["x [N/s^2]","y [N/s^2]","z [N/s^2]"]
    for agent in agent_list:
        u_hist = np.array(agent.history_u)
        t = np.array(agent.history_time)

        if len(u_hist)==0:
            continue
        #x
        axs[0].plot(t,u_hist[:,0],label = agent.agent_id)
        #y
        axs[1].plot(t,u_hist[:,1],label = agent.agent_id)
        #z
        axs[2].plot(t,u_hist[:,2],label = agent.agent_id)

    for i in range(3):
        axs[i].set_ylabel(titles[i])
        axs[i].grid(True)
        axs[i].legend()

    axs[2].set_xlabel("Time [s]")
    fig.suptitle("UAV Control Input u Over Time", fontsize=14)
    plt.tight_layout()

    # === 保存图片 ===
    fig.savefig(os.path.join(output_dir,"u.pdf"), dpi=300, bbox_inches="tight")

def plot_u_formation_control(agent_list):
    fig,axs = plt.subplots(3,1,figsize = (10,10),sharex = True)

    titles = ["x [N/s^2]","y [N/s^2]","z [N/s^2]"]
    for agent in agent_list:
        u_formation_control_hist = np.array(agent.history_u_formation_control)
        
        t = np.array(agent.history_time)

        if len(u_formation_control_hist)==0:
            continue
        #x
        axs[0].plot(t,u_formation_control_hist[:,0],label = agent.agent_id)
        #y
        axs[1].plot(t,u_formation_control_hist[:,1],label = agent.agent_id)
        #z
        axs[2].plot(t,u_formation_control_hist[:,2],label = agent.agent_id)

    for i in range(3):
        axs[i].set_ylabel(titles[i])
        axs[i].grid(True)
        axs[i].legend()

    axs[2].set_xlabel("Time [s]")
    fig.suptitle("UAV Formation Control Input Over Time", fontsize=14)
    plt.tight_layout()

    # === 保存图片 ===
    fig.savefig(os.path.join(output_dir,"u_formation_control.pdf"), dpi=300, bbox_inches="tight")


def plot_gain_evolution(gain_history, dt):
    """
    绘制每个智能体的增益随时间变化曲线
    gain_history: dict 包含每个智能体的增益历史
    dt: float 时间步长
    """
    # 生成时间轴
    # Check if agent_0 exists
    if "agent_0" not in gain_history:
        return
        
    agent_0_angle_gain = gain_history["agent_0"]["angle_gain"]
    N = len(agent_0_angle_gain)
    t = np.arange(N) * dt
    
    # 绘制agent 0的增益曲线（包含距离增益）
    fig, axs = plt.subplots(5, 1, figsize=(10, 12), sharex=True)
    fig.suptitle("Agent 0 Gain Evolution Over Time", fontsize=14)
    
    # 角度增益
    axs[0].plot(t, gain_history["agent_0"]["angle_gain"], label="angle_gain")
    axs[0].set_ylabel("Angle Gain")
    axs[0].grid(True)
    axs[0].legend()
    
    # 距离增益
    if "distance_gain_12" in gain_history["agent_0"]:
        axs[1].plot(t, gain_history["agent_0"]["distance_gain_12"], label="distance_gain_12")
    if "distance_gain_14" in gain_history["agent_0"]:
        axs[1].plot(t, gain_history["agent_0"]["distance_gain_14"], label="distance_gain_14")
    axs[1].set_ylabel("Distance Gains")
    axs[1].grid(True)
    axs[1].legend()
    
    # 阻尼增益x
    axs[2].plot(t, gain_history["agent_0"]["damp_gain_x"], label="damp_gain_x")
    axs[2].set_ylabel("Damp Gain X")
    axs[2].grid(True)
    axs[2].legend()
    
    # 阻尼增益y
    axs[3].plot(t, gain_history["agent_0"]["damp_gain_y"], label="damp_gain_y")
    axs[3].set_ylabel("Damp Gain Y")
    axs[3].grid(True)
    axs[3].legend()
    
    axs[4].plot(t, gain_history["agent_0"]["damp_gain_z"], label="damp_gain_z")
    axs[4].set_ylabel("Damp Gain Z")
    axs[4].grid(True)
    axs[4].legend()

    axs[4].set_xlabel("Time [s]")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,"agent_0_gain_evolution.pdf"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    
    # 绘制agent 1-3的增益曲线（不包含距离增益）
    for agent_id in range(1, 4):
        key = f"agent_{agent_id}"
        if key not in gain_history:
            continue
            
        fig, axs = plt.subplots(4, 1, figsize=(10, 9), sharex=True)
        fig.suptitle(f"Agent {agent_id} Gain Evolution Over Time", fontsize=14)
        
        # 角度增益
        axs[0].plot(t, gain_history[key]["angle_gain"], label="angle_gain")
        axs[0].set_ylabel("Angle Gain")
        axs[0].grid(True)
        axs[0].legend()
        
        # 阻尼增益x
        axs[1].plot(t, gain_history[key]["damp_gain_x"], label="damp_gain_x")
        axs[1].set_ylabel("Damp Gain X")
        axs[1].grid(True)
        axs[1].legend()
        
        # 阻尼增益y
        axs[2].plot(t, gain_history[key]["damp_gain_y"], label="damp_gain_y")
        axs[2].set_ylabel("Damp Gain Y")
        axs[2].grid(True)
        axs[2].legend()
        
        axs[3].plot(t, gain_history[key]["damp_gain_z"], label="damp_gain_z")
        axs[3].set_ylabel("Damp Gain Z")
        axs[3].grid(True)
        axs[3].legend()

        axs[3].set_xlabel("Time [s]")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir,f"agent_{agent_id}_gain_evolution.pdf"), dpi=300, bbox_inches="tight")
        plt.close(fig)
    
def plot_vol_error_evolution(vol_error_list, dt):
    """
    绘制有向四面体体积误差随时间的变化曲线
    :param vol_error_list: list[float], 体积误差历史数据
    :param dt: float, 采样时间步长
    """
    if not vol_error_list:
        print("Warning: vol_error_list is empty, skipping plot.")
        return

    # 生成时间轴
    t = np.arange(len(vol_error_list)) * dt

    plt.figure(figsize=(10, 6))
    
    # 绘制误差曲线
    plt.plot(t, vol_error_list, label="Signed Volume ($V_{1234}$)", color='purple', linewidth=1.5)
    
    # 添加 y=0 的参考线，方便观察收敛情况
    plt.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.6, label="Zero Reference")

    plt.xlabel("Time [s]")
    plt.ylabel("Volume Error [$m^3$]")
    plt.title("Signed Tetrahedron Volume Error Evolution")
    plt.legend()
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    save_path = os.path.join(output_dir, "volume_error_evolution.pdf")
    plt.savefig(save_path, dpi=300, bbox_inches="tight")

def plot_norm_error_evolution(norm_error_list, dt):
    """
    绘制有向四面体体积误差随时间的变化曲线
    :param vol_error_list: list[float], 体积误差历史数据
    :param dt: float, 采样时间步长
    """
    if not norm_error_list:
        print("Warning: norm_error_list is empty, skipping plot.")
        return

    # 生成时间轴
    t = np.arange(len(norm_error_list)) * dt

    plt.figure(figsize=(10, 6))
    
    # 绘制误差曲线
    plt.plot(t, norm_error_list, label="norm vector error()", color='purple', linewidth=1.5)
    
    # 添加 y=0 的参考线，方便观察收敛情况
    plt.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.6, label="Zero Reference")

    plt.xlabel("Time [s]")
    plt.ylabel("Normal Vector Error")
    plt.title("Normal Vector Error Evolution")
    plt.legend()
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    save_path = os.path.join(output_dir, "norm_error_evolution.pdf")
    plt.savefig(save_path, dpi=300, bbox_inches="tight")

def plot_disturbance(agent_list, dt):
    """
    绘制预测扰动和实际扰动的变化曲线
    :param agent_list: 代理列表
    :param dt: 时间步长 (0.002s)
    """
    control_dt = 1/50  # 控制频率50Hz
    num_control_steps = len(agent_list[0].controller.history_dist_hat_acc)
    if num_control_steps == 0:
        return

    control_indices = [0]
    if num_control_steps > 1:
        control_indices.extend(range(9, 9 + 10 * (num_control_steps - 1), 10))
    time = np.array([dt] + [k * control_dt for k in range(1, num_control_steps)])
    
    for idx, agent in enumerate(agent_list):
        dist_acc_actual_full = np.array(agent.history_dist_acc_actual)
        dist_acc_hat = np.array(agent.controller.history_dist_hat_acc)
        
        sampled_indices = [i for i in control_indices if i < len(dist_acc_actual_full)]
        dist_acc_actual = dist_acc_actual_full[sampled_indices][:num_control_steps]
        dist_acc_hat = dist_acc_hat[:len(dist_acc_actual)]
        time_axis = time[:len(dist_acc_actual)]
        
        fig, axes = plt.subplots(3, 1, figsize=(10, 8))
        axes = axes.flatten()
        # NED → ENU: East=NED_Y[1], North=NED_X[0], Up=-NED_Z[2]
        axis_labels = ['X', 'Y', 'Z']
        ned_idx = [1, 0, 2]
        ned_sign = [1, 1, -1]

        for i in range(3):
            val_act = ned_sign[i] * dist_acc_actual[:, ned_idx[i]]
            val_hat = ned_sign[i] * dist_acc_hat[:, ned_idx[i]]
            axes[i].plot(time_axis, val_act, label='Actual Disturbance', color='red')
            observer_label = getattr(agent.controller, "observer_mode_label", "Observer")
            axes[i].plot(
                time_axis,
                val_hat,
                label=f'Estimated Disturbance ({observer_label})',
                color='blue',
                linestyle='--',
            )
            axes[i].set_title(f'Agent {idx+1} - {axis_labels[i]} Disturbance')
            axes[i].set_xlabel('Time (s)')
            axes[i].set_ylabel('Disturbance (m/s²)')
            axes[i].legend()
            axes[i].grid(True)
        
        plt.tight_layout()
        save_path = os.path.join(output_dir, f'disturbance_agent_{idx+1}.png')
        plt.savefig(save_path, dpi=300, bbox_inches="tight")


def plot_disturbance_by_axis(agent_list, dt):
    """
    Generate 3 figures (East, North, Up axes in ENU frame), each showing
    LESO disturbance estimation for all 4 UAVs in a 2x2 subplot layout.
    Disturbance data is converted from NED to ENU to match trajectory plots.
    """
    control_dt = 1/50
    num_control_steps = len(agent_list[0].controller.history_dist_hat_acc)
    if num_control_steps == 0:
        return

    control_indices = [0]
    if num_control_steps > 1:
        control_indices.extend(range(9, 9 + 10 * (num_control_steps - 1), 10))
    time = np.array([dt] + [k * control_dt for k in range(1, num_control_steps)])

    # NED → ENU: East=NED_Y[1], North=NED_X[0], Up=-NED_Z[2]
    axis_labels = ['X', 'Y', 'Z']
    ned_idx = [1, 0, 2]
    ned_sign = [1, 1, -1]

    for axis in range(3):
        fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.5))
        axes = axes.flatten()

        for idx, agent in enumerate(agent_list):
            dist_acc_actual_full = np.array(agent.history_dist_acc_actual)
            dist_acc_hat = np.array(agent.controller.history_dist_hat_acc)

            sampled_indices = [i for i in control_indices if i < len(dist_acc_actual_full)]
            dist_acc_actual = dist_acc_actual_full[sampled_indices][:num_control_steps]
            dist_acc_hat = dist_acc_hat[:len(dist_acc_actual)]
            time_axis = time[:len(dist_acc_actual)]

            ax = axes[idx]
            val_act = ned_sign[axis] * dist_acc_actual[:, ned_idx[axis]]
            val_hat = ned_sign[axis] * dist_acc_hat[:, ned_idx[axis]]
            ax.plot(time_axis, val_act,
                    label='Actual', color='red', linewidth=0.8)
            ax.plot(time_axis, val_hat,
                    label='Estimated', color='blue', linestyle='--', linewidth=0.8)
            ax.set_title(f'QUAV \\#{idx+1}')
            ax.set_xlabel('Time (s)')
            ax.set_ylabel(f'{axis_labels[axis]} (m/s$^2$)')
            ax.legend(fontsize='small')
            ax.grid(True)

        plt.tight_layout()
        save_path = os.path.join(output_dir,
                                 f'fig_exp2_leso_{axis_labels[axis].lower()}.pdf')
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()



