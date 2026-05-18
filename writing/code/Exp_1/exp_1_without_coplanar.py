from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def unit(vec: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    norm = np.linalg.norm(vec)
    if norm < eps:
        return np.zeros_like(vec)
    return vec / norm


def saturation_gain(error: float, threshold: float, bp: float, mu: float) -> float:
    mag = abs(error)
    if mag > threshold:
        return bp * (mag ** (mu - 1.0))
    return bp * (threshold ** (mu - 1.0))


def angle_at_i(
    pi: np.ndarray,
    pj: np.ndarray,
    pk: np.ndarray,
    eps: float = 1e-12,
) -> float:
    rij = pj - pi
    rik = pk - pi
    nij = np.linalg.norm(rij)
    nik = np.linalg.norm(rik)
    if nij < eps or nik < eps:
        return 0.0

    rij = rij / nij
    rik = rik / nik
    dot = float(np.clip(np.dot(rij, rik), -1.0, 1.0))
    cross_norm = float(np.linalg.norm(np.cross(rij, rik)))
    return float(math.atan2(cross_norm, dot))


# ═══════════════════════════════════════════════════════════════════════
#  Desired formation parameters (matching paper Section V-B)
# ═══════════════════════════════════════════════════════════════════════
#  ℓ* = 1.0 m (side length),  n_d = [0, 0, -1]^T (horizontal plane)
#  v_c* = [0.1, 0.1, 0]^T m/s,  α* = π/2 (four right angles)

def get_v_star(t: float) -> np.ndarray:
    return np.array([0.1, 0.1, 0.0], dtype=float)


def get_desired_distance(t: float) -> float:
    return 1.0


def get_desired_normal(v_star: np.ndarray, t: float) -> np.ndarray:
    return np.array([0.0, 0.0, -1.0], dtype=float)


@dataclass
class SimulationConfig:
    """Configuration for Experiment 1 — reduced controller (k_cop = 0).

    Effective proportional gains near equilibrium (paper Sec. V-B):
      k_θ ≈ 1.67,  k_d ≈ 3.70,  k_n ≈ 4.01,  k_v ≈ 5.66
    """
    sim_dt: float = 1.0 / 500.0
    control_dt: float = 1.0 / 50.0
    total_time: float = 150.0
    use_saturation_control: bool = True
    save_animation: bool = True
    gif_frame_stride: int = 120
    gif_fps: int = 40
    undirected_normal_hysteresis: float = 0.05
    output_dir: Path = ROOT / "output" / "double_integrator" / "exp_1_without_coplanar"


class DoubleIntegratorAgent:
    def __init__(self, agent_id: int, p0: list[float], v0: list[float] | None = None):
        self.id = agent_id
        self.pos = np.asarray(p0, dtype=float).reshape(3)
        self.vel = np.zeros(3, dtype=float) if v0 is None else np.asarray(v0, dtype=float).reshape(3)

        self.last_formation_cmd = np.zeros(3, dtype=float)
        self.last_total_acc = np.zeros(3, dtype=float)

        self.pos_hist = [self.pos.copy()]
        self.vel_hist = [self.vel.copy()]
        self.formation_cmd_hist: list[np.ndarray] = []
        self.total_acc_hist: list[np.ndarray] = []

    def step(
        self,
        dt: float,
        formation_cmd: np.ndarray,
    ) -> None:
        total_acc = np.asarray(formation_cmd, dtype=float).reshape(3)

        self.pos = self.pos + self.vel * dt + 0.5 * total_acc * dt * dt
        self.vel = self.vel + total_acc * dt

        self.last_formation_cmd = total_acc.copy()
        self.last_total_acc = total_acc.copy()

        self.pos_hist.append(self.pos.copy())
        self.vel_hist.append(self.vel.copy())
        self.formation_cmd_hist.append(self.last_formation_cmd.copy())
        self.total_acc_hist.append(self.last_total_acc.copy())


class UndirectedFormationController:
    """Reduced four-term controller (paper eq. 40 with k_cop = 0).

    Gains use nonlinear saturation (see Sec. V): k(e) = bp * d_p^(mu-1)
    for |e| > d_p, constant k = bp * d_p^(mu-1) otherwise.
    """

    def __init__(self, config: SimulationConfig):
        self.config = config
        self.use_saturation_control = bool(config.use_saturation_control)
        self.undirected_normal_hysteresis = max(0.0, float(config.undirected_normal_hysteresis))
        self.normal_branch_sign = 1.0

        # ── saturation gain parameters ─────────────────────────
        #  equilibrium effective gains (paper Sec. V-B):
        #    k_θ ≈ 1.67, k_d ≈ 3.70, k_n ≈ 4.01, k_v ≈ 5.66

        # angle term: k_θ(e) · e · (z_ij + z_ik)
        self.bp_angle = 0.5
        self.d_p_angle = 0.35
        self.mu_p_angle = 0.3

        # distance term: k_d(e) · e · z_ij
        self.bp_dist = 0.3
        self.d_p_dist = 0.2
        self.mu_p_dist = 0.3

        # velocity damping: k_v(e) · e
        self.bp_damp = 1.0
        self.d_d_damp = 0.15
        self.mu_d_damp = 0.15

        # normal (attitude) term: k_n(e) · (e_n × r_i)
        self.bp_att = 0.4
        self.d_p_att = 0.1
        self.mu_p_att = 0.3

        # coplanarity term is INTENTIONALLY DISABLED for ablation

        self.angle_constraints = [
            {"name": "e_1", "i": 0, "j": 1, "k": 3, "theta_star": math.pi / 2.0},
            {"name": "e_2", "i": 1, "j": 2, "k": 0, "theta_star": math.pi / 2.0},
            {"name": "e_3", "i": 2, "j": 3, "k": 1, "theta_star": math.pi / 2.0},
            {"name": "e_4", "i": 3, "j": 0, "k": 2, "theta_star": math.pi / 2.0},
        ]
        self.edge_constraints = [
            {"name": "edge12", "i": 0, "j": 1},
            {"name": "edge23", "i": 1, "j": 2},
            {"name": "edge34", "i": 2, "j": 3},
            {"name": "edge41", "i": 3, "j": 0},
        ]

        self.control_time_history: list[float] = []
        self.angle_error_history = {constraint["name"]: [] for constraint in self.angle_constraints}
        self.edge_distance_history = {constraint["name"]: [] for constraint in self.edge_constraints}
        self.coplanar_error_history: list[float] = []
        self.normal_error_history: list[float] = []
        self.normal_selected_alignment_history: list[float] = []

    def _select_undirected_normal(
        self,
        current_normal: np.ndarray,
        desired_normal: np.ndarray,
    ) -> np.ndarray:
        alignment = float(np.dot(current_normal, desired_normal))
        hysteresis = self.undirected_normal_hysteresis
        if alignment > hysteresis:
            self.normal_branch_sign = 1.0
        elif alignment < -hysteresis:
            self.normal_branch_sign = -1.0
        return self.normal_branch_sign * desired_normal

    def _gain(
        self,
        error: float,
        threshold: float,
        bp: float,
        mu: float,
    ) -> float:
        """Compute saturation gain (paper eq. 41).

        When use_saturation_control=False, returns constant bp
        (constant-gain mode matching Li et al. (2025) style).
        """
        if self.use_saturation_control:
            return saturation_gain(error, threshold, bp, mu)
        return bp  # constant gain (Li et al. style)

    def control_law_output(
        self,
        agents: list[DoubleIntegratorAgent],
        v_star_stack: np.ndarray,
        t: float,
    ) -> np.ndarray:
        """Evaluate reduced control law (paper eq. 40, k_cop = 0)."""
        u = np.zeros((4, 3), dtype=float)
        p_stack = np.vstack([agent.pos for agent in agents])
        v_stack = np.vstack([agent.vel for agent in agents])

        self.control_time_history.append(t)

        # ── Term 1: angle constraints (paper eq. 40, k_θ term) ──
        for angle in self.angle_constraints:
            i = angle["i"]
            j = angle["j"]
            k = angle["k"]

            zij = unit(p_stack[j] - p_stack[i])
            zik = unit(p_stack[k] - p_stack[i])
            theta = angle_at_i(p_stack[i], p_stack[j], p_stack[k])
            e_angle = theta - angle["theta_star"]

            self.angle_error_history[angle["name"]].append(e_angle)
            gain = self._gain(e_angle, self.d_p_angle, self.bp_angle, self.mu_p_angle)
            u[i] += -gain * e_angle * (zij + zik)

        # ── Log all edge distances ──
        for edge in self.edge_constraints:
            d = float(np.linalg.norm(p_stack[edge["j"]] - p_stack[edge["i"]]))
            self.edge_distance_history[edge["name"]].append(d)

        # ── Term 2: distance constraints (paper eq. 40, k_d term) ──
        desired_edge = get_desired_distance(t)
        d12 = np.linalg.norm(p_stack[1] - p_stack[0])
        d14 = np.linalg.norm(p_stack[3] - p_stack[0])
        z12 = unit(p_stack[1] - p_stack[0])
        z14 = unit(p_stack[3] - p_stack[0])

        e_d12 = d12 - desired_edge
        e_d14 = d14 - desired_edge
        gain_d12 = self._gain(e_d12, self.d_p_dist, self.bp_dist, self.mu_p_dist)
        gain_d14 = self._gain(e_d14, self.d_p_dist, self.bp_dist, self.mu_p_dist)

        u[0] += gain_d12 * e_d12 * z12
        u[0] += gain_d14 * e_d14 * z14
        u[1] -= gain_d12 * e_d12 * z12
        u[3] -= gain_d14 * e_d14 * z14

        # ── Term 5: velocity damping (paper eq. 40, k_v term) ──
        e_v_stack = v_stack - v_star_stack
        for i in range(4):
            for axis in range(3):
                e_v = e_v_stack[i, axis]
                gain = self._gain(e_v, self.d_d_damp, self.bp_damp, self.mu_d_damp)
                u[i, axis] += -gain * e_v

        # ── Coplanarity term DISABLED for ablation (k_cop = 0) ──
        z31 = unit(p_stack[0] - p_stack[2])
        z32 = unit(p_stack[1] - p_stack[2])
        z34 = unit(p_stack[3] - p_stack[2])
        n234_local = unit(np.cross(z32, z34))
        coplanar_error = float(np.dot(z31, n234_local))
        self.coplanar_error_history.append(coplanar_error)
        # coplanarity force deliberately omitted — see paper Sec. V-A

        # ── Term 3: undirected plane-normal control (paper eq. 40, k_n term) ──
        n142 = unit(np.cross(z14, z12))                     # n_{124} = (b12×b14)/‖·‖
        desired_normal = unit(get_desired_normal(v_star_stack[0], t))
        selected_normal = unit(self._select_undirected_normal(n142, desired_normal))
        normal_error = np.cross(n142, selected_normal)      # e_n = n_{124} × ̄n_d
        undirected_angle = math.acos(float(np.clip(abs(np.dot(n142, desired_normal)), 0.0, 1.0)))
        self.normal_error_history.append(undirected_angle)
        self.normal_selected_alignment_history.append(float(np.dot(n142, selected_normal)))

        p_center = np.mean(np.vstack([p_stack[0], p_stack[1], p_stack[3]]), axis=0)
        normal_gain = self._gain(
            float(np.linalg.norm(normal_error)),
            self.d_p_att,
            self.bp_att,
            self.mu_p_att,
        )
        for idx in (0, 1, 3):
            r_i = p_stack[idx] - p_center
            u[idx] += normal_gain * np.cross(normal_error, r_i)

        return u


def simulate(config: SimulationConfig) -> tuple[list[DoubleIntegratorAgent], UndirectedFormationController, np.ndarray]:
    """Run Experiment 1 simulation (k_cop = 0, reduced controller).

    Initial conditions (paper eq. 43):
      p1(0)=[0, 0.6, 0]^T,  p2(0)=[0, 0, 0]^T,
      p3(0)=[0.6, 0, 0.05]^T, p4(0)=[0.6, 0.6, 0]^T,
      v_i(0)=[0, 0, 0]^T m/s.
    Agent #3 has a deliberate 0.05 m out-of-plane perturbation.
    """
    agents = [
        DoubleIntegratorAgent(1, [0.0, 0.6, 0.0]),
        DoubleIntegratorAgent(2, [0.0, 0.0, 0.0]),
        DoubleIntegratorAgent(3, [0.6, 0.0, 0.05]),
        DoubleIntegratorAgent(4, [0.6, 0.6, 0.0]),
    ]

    controller = UndirectedFormationController(config)

    t = 0.0
    sim_steps = int(round(config.total_time / config.sim_dt))
    sim_time_hist = np.zeros(sim_steps + 1, dtype=float)
    control_elapsed = 0.0

    v_star_stack = np.tile(get_v_star(t), (4, 1))
    latest_u_formation = controller.control_law_output(agents, v_star_stack, t)

    for step in range(sim_steps):
        for idx, agent in enumerate(agents):
            agent.step(config.sim_dt, latest_u_formation[idx])

        t += config.sim_dt
        sim_time_hist[step + 1] = t
        control_elapsed += config.sim_dt

        if control_elapsed >= config.control_dt - 1e-12:
            while control_elapsed >= config.control_dt - 1e-12:
                control_elapsed -= config.control_dt
            v_star_stack = np.tile(get_v_star(t), (4, 1))
            latest_u_formation = controller.control_law_output(agents, v_star_stack, t)

    return agents, controller, sim_time_hist




def _plot_xy_trajectory(agents: list[DoubleIntegratorAgent], output_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 6))
    for agent in agents:
        pos_hist = np.array(agent.pos_hist)
        ax.plot(pos_hist[:, 0], pos_hist[:, 1], label=f"agent {agent.id}")
        ax.scatter(*pos_hist[-1, :2], color="black", s=20, zorder=5)
        ax.text(
            pos_hist[-1, 0] + 0.04,
            pos_hist[-1, 1] + 0.04,
            str(agent.id),
            fontsize=10,
        )
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title("XY Trajectory")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True)
    ax.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "xy_trajectory.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def _plot_3d_trajectory(agents: list[DoubleIntegratorAgent], output_dir: Path) -> None:
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")
    for agent in agents:
        pos_hist = np.array(agent.pos_hist)
        ax.plot(pos_hist[:, 0], pos_hist[:, 1], pos_hist[:, 2], label=f"agent {agent.id}")
        ax.scatter(*pos_hist[-1], s=30)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_zlabel("z (m)")
    ax.set_title("3D Trajectory")
    ax.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "trajectory_3d.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def _plot_height_over_time(
    agents: list[DoubleIntegratorAgent],
    sim_time: np.ndarray,
    output_dir: Path,
) -> None:
    plt.figure(figsize=(8, 5))
    for agent in agents:
        pos_hist = np.array(agent.pos_hist)
        height = -pos_hist[:, 2]
        plt.plot(sim_time, height, label=f"agent {agent.id}")
    plt.xlabel("time (s)")
    plt.ylabel("height (m)")
    plt.title("Height Over Time")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "z_height.png", dpi=200, bbox_inches="tight")
    plt.close()


def _plot_control_history(
    agents: list[DoubleIntegratorAgent],
    sim_time: np.ndarray,
    output_dir: Path,
    attr: str,
    filename: str,
    title: str,
) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    labels = ["x", "y", "z"]
    time_axis = sim_time[1:]
    for agent in agents:
        hist = np.array(getattr(agent, attr))
        if hist.size == 0:
            continue
        for axis in range(3):
            axes[axis].plot(time_axis, hist[:, axis], label=f"agent {agent.id}")
            axes[axis].set_ylabel(labels[axis])
            axes[axis].grid(True)
    axes[0].set_title(title)
    axes[-1].set_xlabel("time (s)")
    axes[0].legend()
    plt.tight_layout()
    plt.savefig(output_dir / filename, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _plot_speed_and_accel(
    agents: list[DoubleIntegratorAgent],
    sim_time: np.ndarray,
    output_dir: Path,
) -> None:
    time_axis = sim_time[1:]

    plt.figure(figsize=(8, 5))
    for agent in agents:
        vel_hist = np.array(agent.vel_hist[1:])
        speed = np.linalg.norm(vel_hist, axis=1)
        plt.plot(time_axis, speed, label=f"agent {agent.id}")
    plt.xlabel("time (s)")
    plt.ylabel("|v| (m/s)")
    plt.title("Speed")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "speed.png", dpi=200, bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(8, 5))
    for agent in agents:
        acc_hist = np.array(agent.total_acc_hist)
        accel = np.linalg.norm(acc_hist, axis=1)
        plt.plot(time_axis, accel, label=f"agent {agent.id}")
    plt.xlabel("time (s)")
    plt.ylabel("|a| (m/s^2)")
    plt.title("Acceleration")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "accel.png", dpi=200, bbox_inches="tight")
    plt.close()


def _plot_angle_errors(controller: UndirectedFormationController, output_dir: Path) -> None:
    t = np.array(controller.control_time_history)
    plt.figure(figsize=(8, 5))
    for name, values in controller.angle_error_history.items():
        plt.plot(t, values, label=name)
    plt.xlabel("time (s)")
    plt.ylabel("error (rad)")
    plt.title("Angle Errors")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "angle_error.png", dpi=200, bbox_inches="tight")
    plt.close()


def _plot_edge_histories(
    controller: UndirectedFormationController,
    sim_time: np.ndarray,
    output_dir: Path,
) -> None:
    t = np.array(controller.control_time_history)
    edge_target = np.array([get_desired_distance(ti) for ti in t])

    plt.figure(figsize=(8, 5))
    for name, values in controller.edge_distance_history.items():
        plt.plot(t, values, label=name)
    plt.plot(t, edge_target, "--", label="edge target")
    plt.xlabel("time (s)")
    plt.ylabel("distance (m)")
    plt.title("Edge Distances")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "edge_distances.png", dpi=200, bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(8, 5))
    for name, values in controller.edge_distance_history.items():
        plt.plot(t, np.array(values) - edge_target, label=name)
    plt.axhline(0.0, color="black", linewidth=1.0)
    plt.xlabel("time (s)")
    plt.ylabel("error (m)")
    plt.title("Edge Distance Errors")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "edge_distance_errors.png", dpi=200, bbox_inches="tight")
    plt.close()


def _plot_scalar_history(
    time_axis: np.ndarray,
    values: np.ndarray,
    output_dir: Path,
    filename: str,
    title: str,
    ylabel: str,
) -> None:
    plt.figure(figsize=(8, 5))
    plt.plot(time_axis, values)
    plt.xlabel("time (s)")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_dir / filename, dpi=200, bbox_inches="tight")
    plt.close()


def _save_animation(
    agents: list[DoubleIntegratorAgent],
    output_dir: Path,
    sim_dt: float,
    frame_stride: int,
    fps: int,
) -> None:
    pos_hist = np.stack([np.array(agent.pos_hist) for agent in agents], axis=1)
    plot_hist = pos_hist.copy()
    plot_hist[:, :, 2] *= -1.0

    frame_stride = max(1, int(frame_stride))
    fps = max(1, int(fps))
    frames = list(range(0, pos_hist.shape[0], frame_stride))
    if frames[-1] != pos_hist.shape[0] - 1:
        frames.append(pos_hist.shape[0] - 1)
    sampled_hist = plot_hist[frames]

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    all_xyz = plot_hist.reshape(-1, 3)
    margin_xy = 0.5
    xmin, xmax = all_xyz[:, 0].min() - margin_xy, all_xyz[:, 0].max() + margin_xy
    ymin, ymax = all_xyz[:, 1].min() - margin_xy, all_xyz[:, 1].max() + margin_xy

    zmin_raw = all_xyz[:, 2].min()
    zmax_raw = all_xyz[:, 2].max()
    z_center = 0.5 * (zmin_raw + zmax_raw)
    z_half_span = max(0.5, 0.5 * (zmax_raw - zmin_raw) + 0.1)
    zmin, zmax = z_center - z_half_span, z_center + z_half_span

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_zlim(zmin, zmax)
    ax.set_xlabel("East [m]")
    ax.set_ylabel("North [m]")
    ax.set_zlabel("Height [m]")
    ax.set_title("Undirected Normal Control Animation")
    ax.view_init(elev=32, azim=-60)

    colors = ["b", "g", "r", "c"]
    traj_lines = []
    pos_dots = []
    labels = []

    for idx, agent in enumerate(agents):
        color = colors[idx % len(colors)]
        traj_line, = ax.plot([], [], [], linewidth=1.5, color=color, label=f"agent {agent.id}")
        pos_dot, = ax.plot([], [], [], marker="o", markersize=8, color=color)
        label = ax.text(0.0, 0.0, 0.0, "", color="black", fontsize=10)
        traj_lines.append(traj_line)
        pos_dots.append(pos_dot)
        labels.append(label)

    ax.legend()

    def update(frame_idx: int):
        pts = sampled_hist[frame_idx]
        actual_frame = frames[frame_idx]

        for idx in range(len(agents)):
            hist = sampled_hist[: frame_idx + 1, idx, :]
            traj_lines[idx].set_data(hist[:, 0], hist[:, 1])
            traj_lines[idx].set_3d_properties(hist[:, 2])

            pos_dots[idx].set_data([pts[idx, 0]], [pts[idx, 1]])
            pos_dots[idx].set_3d_properties([pts[idx, 2]])

            labels[idx].set_position((pts[idx, 0], pts[idx, 1]))
            labels[idx].set_3d_properties(pts[idx, 2] + 0.04)
            labels[idx].set_text(f"{idx + 1}")

        ax.set_title(f"Undirected Normal Control Animation\nt = {actual_frame * sim_dt:.2f} s")
        return traj_lines + pos_dots + labels

    ani = FuncAnimation(fig, update, frames=len(frames), interval=1000 / fps, blit=False)
    ani.save(output_dir / "formation.gif", writer=PillowWriter(fps=fps))
    plt.close(fig)


def save_outputs(
    agents: list[DoubleIntegratorAgent],
    controller: UndirectedFormationController,
    sim_time: np.ndarray,
    config: SimulationConfig,
) -> None:
    output_dir = config.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    _plot_xy_trajectory(agents, output_dir)
    _plot_3d_trajectory(agents, output_dir)
    _plot_height_over_time(agents, sim_time, output_dir)
    _plot_control_history(
        agents,
        sim_time,
        output_dir,
        "formation_cmd_hist",
        "u_formation_control.png",
        "Formation Control Acceleration",
    )
    _plot_control_history(
        agents,
        sim_time,
        output_dir,
        "total_acc_hist",
        "u_total.png",
        "Total Acceleration",
    )
    _plot_speed_and_accel(agents, sim_time, output_dir)
    _plot_angle_errors(controller, output_dir)
    _plot_edge_histories(controller, sim_time, output_dir)
    _plot_scalar_history(
        np.array(controller.control_time_history),
        np.array(controller.coplanar_error_history),
        output_dir,
        "coplanar_error.png",
        "Coplanar Error",
        "error",
    )
    _plot_scalar_history(
        np.array(controller.control_time_history),
        np.array(controller.normal_error_history),
        output_dir,
        "normal_error.png",
        "Undirected Plane Error",
        "error (rad)",
    )
    _plot_scalar_history(
        np.array(controller.control_time_history),
        np.array(controller.normal_selected_alignment_history),
        output_dir,
        "normal_selected_alignment.png",
        "Selected Normal Alignment",
        r"$n^\top \bar n_d$",
    )
    if config.save_animation:
        _save_animation(
            agents,
            output_dir,
            config.sim_dt,
            frame_stride=config.gif_frame_stride,
            fps=config.gif_fps,
        )


def print_summary(
    agents: list[DoubleIntegratorAgent],
    controller: UndirectedFormationController,
    config: SimulationConfig,
) -> None:
    final_positions = np.vstack([agent.pos for agent in agents])
    final_velocities = np.vstack([agent.vel for agent in agents])

    d12 = np.linalg.norm(final_positions[1] - final_positions[0])
    d14 = np.linalg.norm(final_positions[3] - final_positions[0])
    edge_target = get_desired_distance(config.total_time)

    print(f"Results saved to: {config.output_dir}")
    print("Final positions:")
    print(np.array2string(final_positions, precision=3, suppress_small=True))
    print("Final velocities:")
    print(np.array2string(final_velocities, precision=3, suppress_small=True))
    print(f"Final d12 error: {d12 - edge_target:.4f} m")
    print(f"Final d14 error: {d14 - edge_target:.4f} m")
    print("Final angle errors (rad):")
    final_angle_errors = {
        name: values[-1]
        for name, values in controller.angle_error_history.items()
        if values
    }
    print(final_angle_errors)
    if controller.coplanar_error_history:
        print(f"Final coplanar error: {controller.coplanar_error_history[-1]:.4f}")
    if controller.normal_error_history:
        print(f"Final undirected plane error: {controller.normal_error_history[-1]:.4f}")


def main() -> None:
    config = SimulationConfig()
    agents, controller, sim_time = simulate(config)
    save_outputs(agents, controller, sim_time, config)
    print_summary(agents, controller, config)


if __name__ == "__main__":
    main()
