from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib.animation import FuncAnimation, PillowWriter

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from agent.observer import (  # noqa: E402
    LESO,
)


def _load_uav_square_half_size(default: float = 0.16) -> float:
    params_path = ROOT / "uav_params.yaml"
    try:
        with params_path.open("r", encoding="utf-8") as file:
            params = yaml.safe_load(file)
        arm_length = float(params["uav"]["arm_length"])
    except Exception:
        return default
    return arm_length / math.sqrt(2.0)


UAV_SQUARE_HALF_SIZE = _load_uav_square_half_size()


def unit(vec: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    norm = np.linalg.norm(vec)
    if norm < eps:
        return np.zeros_like(vec)
    return vec / norm


def get_capture_point(
    p_stack: np.ndarray,
    offset_ratio_u: float,
    offset_ratio_v: float,
) -> np.ndarray:
    p_stack = np.asarray(p_stack, dtype=float).reshape(-1, 3)
    net_center = np.mean(p_stack, axis=0)

    left_mid = 0.5 * (p_stack[0] + p_stack[1])
    right_mid = 0.5 * (p_stack[2] + p_stack[3])
    lower_mid = 0.5 * (p_stack[1] + p_stack[2])
    upper_mid = 0.5 * (p_stack[0] + p_stack[3])

    axis_u_vec = right_mid - left_mid
    axis_v_vec = upper_mid - lower_mid
    half_span_u = 0.5 * np.linalg.norm(axis_u_vec)
    half_span_v = 0.5 * np.linalg.norm(axis_v_vec)

    axis_u = unit(axis_u_vec)
    axis_v = unit(axis_v_vec)
    return (
        net_center
        + offset_ratio_u * half_span_u * axis_u
        + offset_ratio_v * half_span_v * axis_v
    )


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


def get_v_star(t: float) -> np.ndarray:
    if t < 30.0:
        return np.array([0.1, 0.1, 0.0], dtype=float)
    if t < 60.0:
        # return np.array([0.1, 0.1, 0.0], dtype=float)
        return np.array([0.2, 0.1, 0.0], dtype=float)
    if t < 100.0:
        return np.array([0.1, 0.1, 0.0], dtype=float)
    return np.array([0.1, 0.1, 0.0], dtype=float)


def get_desired_distance(t: float) -> float:
    # return 1.0 if t < 100.0 else 0.64
    # if 90<=t<=100:
    #     return 0.34
    return 1.0


def get_desired_normal(v_star: np.ndarray, t: float) -> np.ndarray:
    v_xy = np.array([v_star[0], v_star[1], 0.0], dtype=float)
    v_norm = np.linalg.norm(v_xy)
    if v_norm < 1e-3:
        return np.array([0.0, 0.0, -1.0], dtype=float)

    n_forward = v_xy / v_norm
    # if t >= 60.0:
    #     return unit(-n_forward + np.array([0.0, 0.0, -1.0]))
    # if t < 30.0:
    #     return unit(-n_forward + np.array([0.0, 0.0, 1.0]))
    if 30.0 <= t <= 90.0:
        return n_forward
    # if 60.0<= t <=90:
    #     return unit(-n_forward + np.array([0.0, 0.0, 1.0]))
    return unit(np.array([0.0, 0.0, 1.0], dtype=float))


def get_disturbance(
    t: float,
    agent_idx: int,
    p_stack: np.ndarray,
    config: "SimulationConfig",
) -> np.ndarray:
    disturbance = np.zeros(3, dtype=float)
    if not (config.capture_start_time <= t <= config.capture_end_time):
        return disturbance

    p_stack = np.asarray(p_stack, dtype=float).reshape(-1, 3)
    if not 0 <= agent_idx < len(p_stack):
        return disturbance

    # Off-center capture: the target impact point stays inside the current net
    # but is shifted away from the geometric center.
    capture_point = get_capture_point(
        p_stack,
        config.capture_point_offset_ratio[0],
        config.capture_point_offset_ratio[1],
    )
    pull_direction = unit(capture_point - p_stack[agent_idx])
    disturbance = (config.capture_pull_force / config.mass) * pull_direction
    return disturbance


@dataclass
class SimulationConfig:
    sim_dt: float = 1.0 / 500.0
    control_dt: float = 1.0 / 50.0
    total_time: float = 200.0
    mass: float = 1.5
    capture_pull_force: float = 1.0
    capture_start_time: float = 60.0
    capture_end_time: float = 65.0
    capture_point_offset_ratio: tuple[float, float] = (0.25, -0.15)
    use_leso: bool = True
    use_saturation_control: bool = True
    save_animation: bool = True
    gif_frame_stride: int = 120
    gif_fps: int = 40
    undirected_normal_hysteresis: float = 0.05
    output_dir: Path = ROOT / "output" / "double_integrator" / "current_formation_cotrol_undirected"


class DoubleIntegratorAgent:
    def __init__(self, agent_id: int, p0: list[float], v0: list[float] | None = None):
        self.id = agent_id
        self.pos = np.asarray(p0, dtype=float).reshape(3)
        self.vel = np.zeros(3, dtype=float) if v0 is None else np.asarray(v0, dtype=float).reshape(3)

        self.last_formation_cmd = np.zeros(3, dtype=float)
        self.last_total_acc = np.zeros(3, dtype=float)
        self.last_rope_force = np.zeros(3, dtype=float)
        self.last_disturbance = np.zeros(3, dtype=float)

        self.pos_hist = [self.pos.copy()]
        self.vel_hist = [self.vel.copy()]
        self.formation_cmd_hist: list[np.ndarray] = []
        self.total_acc_hist: list[np.ndarray] = []
        self.rope_force_hist: list[np.ndarray] = []
        self.disturbance_hist: list[np.ndarray] = []

    def step(
        self,
        pos_next: np.ndarray,
        vel_next: np.ndarray,
        formation_cmd: np.ndarray,
        rope_force: np.ndarray,
        disturbance: np.ndarray,
        mass: float,
    ) -> None:
        pos_next = np.asarray(pos_next, dtype=float).reshape(3)
        vel_next = np.asarray(vel_next, dtype=float).reshape(3)
        formation_cmd = np.asarray(formation_cmd, dtype=float).reshape(3)
        rope_force = np.asarray(rope_force, dtype=float).reshape(3)
        disturbance = np.asarray(disturbance, dtype=float).reshape(3)

        total_acc = formation_cmd + rope_force / mass + disturbance

        self.pos = pos_next
        self.vel = vel_next

        self.last_formation_cmd = formation_cmd.copy()
        self.last_total_acc = total_acc.copy()
        self.last_rope_force = rope_force.copy()
        self.last_disturbance = disturbance.copy()

        self.pos_hist.append(self.pos.copy())
        self.vel_hist.append(self.vel.copy())
        self.formation_cmd_hist.append(self.last_formation_cmd.copy())
        self.total_acc_hist.append(self.last_total_acc.copy())
        self.rope_force_hist.append(self.last_rope_force.copy())
        self.disturbance_hist.append(self.last_disturbance.copy())


class UndirectedFormationController:
    def __init__(self, config: SimulationConfig):
        self.config = config
        self.use_leso = config.use_leso
        self.use_saturation_control = bool(config.use_saturation_control)
        self.undirected_normal_hysteresis = max(0.0, float(config.undirected_normal_hysteresis))
        self.normal_branch_sign = 1.0

        self.K_rope = 74.0
        self.C_rope = 18.0
        self.l0 = 1.0

        self.bp_angle = 1.2#0.8
        self.d_p_angle = 0.35
        self.mu_p_angle = 0.3

        self.bp_dist = 1.2
        self.d_p_dist = 0.2
        self.mu_p_dist = 0.3

        self.bp_damp = 1.3
        self.d_d_damp = 0.15
        self.mu_d_damp = 0.15

        self.bp_att = 0.8
        self.d_p_att = 0.1
        self.mu_p_att = 0.3

        self.bp_cop = 1.0
        self.d_p_cop = 0.05
        self.mu_p_cop = 0.5

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
        self.tension_history = {constraint["name"]: [] for constraint in self.edge_constraints}
        self.coplanar_error_history: list[float] = []
        self.normal_error_history: list[float] = []
        self.normal_selected_alignment_history: list[float] = []

        self.leso_dt = config.control_dt
        self.leso_omega = 22 #15.0

        self.observers: list[dict[str, object]] = []
        self.last_applied_acc_cmds: list[np.ndarray | None] = []
        self.dist_hat_history = {f"agent_{i}": [] for i in range(4)}

    def _build_observer_axis(self, initial_p: float, initial_v: float):
        omega_o = self.leso_omega
        return LESO(
            beta1=3.0 * omega_o,
            beta2=3.0 * omega_o**2,
            beta3=omega_o**3,
            dt=self.leso_dt,
            initial_p=initial_p,
            initial_v=initial_v,
        )

    def reset_observers(self, agents: list[DoubleIntegratorAgent]) -> None:
        if not self.use_leso:
            return

        self.observers = []
        self.last_applied_acc_cmds = []
        for agent in agents:
            self.observers.append(
                {
                    "x": self._build_observer_axis(agent.pos[0], agent.vel[0]),
                    "y": self._build_observer_axis(agent.pos[1], agent.vel[1]),
                    "z": self._build_observer_axis(agent.pos[2], agent.vel[2]),
                }
            )
            self.last_applied_acc_cmds.append(None)

    def _compensate_with_leso(
        self,
        agents: list[DoubleIntegratorAgent],
        u_raw: np.ndarray,
    ) -> np.ndarray:
        if not self.use_leso:
            return np.array(u_raw, dtype=float, copy=True)

        if len(self.observers) != len(agents):
            self.reset_observers(agents)

        u_comp = np.array(u_raw, dtype=float, copy=True)
        for idx, agent in enumerate(agents):
            prev_u = self.last_applied_acc_cmds[idx]
            if prev_u is None:
                prev_u = np.array(u_raw[idx], dtype=float, copy=True)

            axes = self.observers[idx]
            _, _, dist_hat_x = axes["x"].update(agent.pos[0], prev_u[0])
            _, _, dist_hat_y = axes["y"].update(agent.pos[1], prev_u[1])
            _, _, dist_hat_z = axes["z"].update(agent.pos[2], prev_u[2])

            dist_hat = np.array([dist_hat_x, dist_hat_y, dist_hat_z], dtype=float)
            self.dist_hat_history[f"agent_{idx}"].append(dist_hat.copy())
            u_comp[idx] -= dist_hat
            self.last_applied_acc_cmds[idx] = u_comp[idx].copy()

        return u_comp

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
        if self.use_saturation_control:
            return saturation_gain(error, threshold, bp, mu)
        return bp

    def control_law_output(
        self,
        agents: list[DoubleIntegratorAgent],
        v_star_stack: np.ndarray,
        t: float,
    ) -> np.ndarray:
        u = np.zeros((4, 3), dtype=float)
        p_stack = np.vstack([agent.pos for agent in agents])
        v_stack = np.vstack([agent.vel for agent in agents])

        self.control_time_history.append(t)

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

        e_v_stack = v_stack - v_star_stack
        for i in range(4):
            for axis in range(3):
                e_v = e_v_stack[i, axis]
                gain = self._gain(e_v, self.d_d_damp, self.bp_damp, self.mu_d_damp)
                u[i, axis] += -gain * e_v

        z31 = unit(p_stack[0] - p_stack[2])
        z32 = unit(p_stack[1] - p_stack[2])
        z34 = unit(p_stack[3] - p_stack[2])
        n234_local = unit(np.cross(z32, z34))
        coplanar_error = float(np.dot(z31, n234_local))
        self.coplanar_error_history.append(coplanar_error)

        cop_gain = self._gain(coplanar_error, self.d_p_cop, self.bp_cop, self.mu_p_cop)
        u[2] -= cop_gain * coplanar_error * n234_local

        n142 = unit(np.cross(z14, z12))
        desired_normal = unit(get_desired_normal(v_star_stack[0], t))
        selected_normal = unit(self._select_undirected_normal(n142, desired_normal))
        normal_error = np.cross(n142, selected_normal)
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

        return self._compensate_with_leso(agents, u)

    def _tension_law_from_state(
        self,
        p_stack: np.ndarray,
        v_stack: np.ndarray,
        *,
        record_history: bool,
    ) -> np.ndarray:
        p_stack = np.asarray(p_stack, dtype=float).reshape(-1, 3)
        v_stack = np.asarray(v_stack, dtype=float).reshape(-1, 3)
        rope_force = np.zeros_like(p_stack)

        for edge in self.edge_constraints:
            i = edge["i"]
            j = edge["j"]
            rij = p_stack[j] - p_stack[i]
            distance = np.linalg.norm(rij)
            if record_history:
                self.edge_distance_history[edge["name"]].append(distance)

            tension = 0.0
            if distance > 1e-12:
                zij = rij / distance
                if distance > self.l0:
                    d_dot = float(np.dot(v_stack[j] - v_stack[i], zij))
                    tension = self.K_rope * (distance - self.l0) + self.C_rope * d_dot
                    tension = max(0.0, tension)
                    rope_force[i] += tension * zij
                    rope_force[j] -= tension * zij

            if record_history:
                self.tension_history[edge["name"]].append(tension)

        return rope_force

    def tension_law_output(self, agents: list[DoubleIntegratorAgent]) -> np.ndarray:
        p_stack = np.vstack([agent.pos for agent in agents])
        v_stack = np.vstack([agent.vel for agent in agents])
        return self._tension_law_from_state(p_stack, v_stack, record_history=True)


def _stack_agent_state(agents: list[DoubleIntegratorAgent]) -> np.ndarray:
    p_stack = np.vstack([agent.pos for agent in agents])
    v_stack = np.vstack([agent.vel for agent in agents])
    return np.concatenate((p_stack.reshape(-1), v_stack.reshape(-1)))


def _split_agent_state(state: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    state = np.asarray(state, dtype=float)
    half = state.size // 2
    p_stack = state[:half].reshape(-1, 3)
    v_stack = state[half:].reshape(-1, 3)
    return p_stack, v_stack


def _disturbance_stack(
    t: float,
    p_stack: np.ndarray,
    config: SimulationConfig,
) -> np.ndarray:
    p_stack = np.asarray(p_stack, dtype=float).reshape(-1, 3)
    return np.vstack(
        [get_disturbance(t, idx, p_stack, config) for idx in range(len(p_stack))]
    )


def _formation_dynamics(
    state: np.ndarray,
    t: float,
    formation_cmd: np.ndarray,
    controller: UndirectedFormationController,
    config: SimulationConfig,
) -> np.ndarray:
    p_stack, v_stack = _split_agent_state(state)
    formation_cmd = np.asarray(formation_cmd, dtype=float).reshape(p_stack.shape)
    rope_force = controller._tension_law_from_state(
        p_stack,
        v_stack,
        record_history=False,
    )
    disturbance = _disturbance_stack(t, p_stack, config)
    acc = formation_cmd + rope_force / config.mass + disturbance
    return np.concatenate((v_stack.reshape(-1), acc.reshape(-1)))


def _rk4_integrate_step(
    state: np.ndarray,
    t: float,
    dt: float,
    formation_cmd: np.ndarray,
    controller: UndirectedFormationController,
    config: SimulationConfig,
) -> np.ndarray:
    k1 = _formation_dynamics(state, t, formation_cmd, controller, config)
    k2 = _formation_dynamics(
        state + 0.5 * dt * k1,
        t + 0.5 * dt,
        formation_cmd,
        controller,
        config,
    )
    k3 = _formation_dynamics(
        state + 0.5 * dt * k2,
        t + 0.5 * dt,
        formation_cmd,
        controller,
        config,
    )
    k4 = _formation_dynamics(
        state + dt * k3,
        t + dt,
        formation_cmd,
        controller,
        config,
    )
    return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def simulate(config: SimulationConfig) -> tuple[list[DoubleIntegratorAgent], UndirectedFormationController, np.ndarray]:
    agents = [
        DoubleIntegratorAgent(1, [0.0, 0.34, 0.5]),
        DoubleIntegratorAgent(2, [0.0, 0.0, -0.2]),
        DoubleIntegratorAgent(3, [0.34, 0.0, 0.1]),
        DoubleIntegratorAgent(4, [0.34, 0.34, 0.3]),
    ]

    controller = UndirectedFormationController(config)
    controller.reset_observers(agents)

    t = 0.0
    sim_steps = int(round(config.total_time / config.sim_dt))
    sim_time_hist = np.zeros(sim_steps + 1, dtype=float)
    control_elapsed = 0.0

    v_star_stack = np.tile(get_v_star(t), (4, 1))
    latest_u_formation = controller.control_law_output(agents, v_star_stack, t)

    for step in range(sim_steps):
        state = _stack_agent_state(agents)
        p_stack, v_stack = _split_agent_state(state)
        rope_force = controller._tension_law_from_state(
            p_stack,
            v_stack,
            record_history=True,
        )
        disturbance_stack = _disturbance_stack(t, p_stack, config)
        next_state = _rk4_integrate_step(
            state,
            t,
            config.sim_dt,
            latest_u_formation,
            controller,
            config,
        )
        next_p_stack, next_v_stack = _split_agent_state(next_state)
        for idx, agent in enumerate(agents):
            agent.step(
                next_p_stack[idx],
                next_v_stack[idx],
                latest_u_formation[idx],
                rope_force[idx],
                disturbance_stack[idx],
                config.mass,
            )

        t += config.sim_dt
        sim_time_hist[step + 1] = t
        control_elapsed += config.sim_dt

        if control_elapsed >= config.control_dt - 1e-12:
            while control_elapsed >= config.control_dt - 1e-12:
                control_elapsed -= config.control_dt
            v_star_stack = np.tile(get_v_star(t), (4, 1))
            latest_u_formation = controller.control_law_output(agents, v_star_stack, t)

    return agents, controller, sim_time_hist


def _square_outline(center_xy: np.ndarray, size: float = 0.03) -> np.ndarray:
    cx, cy = np.asarray(center_xy, dtype=float).reshape(2)
    return np.array(
        [
            [cx + size, cy + size],
            [cx + size, cy - size],
            [cx - size, cy - size],
            [cx - size, cy + size],
            [cx + size, cy + size],
        ],
        dtype=float,
    )


def _square_outline_3d(center_xyz: np.ndarray, size: float = 0.03) -> np.ndarray:
    cx, cy, cz = np.asarray(center_xyz, dtype=float).reshape(3)
    return np.array(
        [
            [cx + size, cy + size, cz],
            [cx + size, cy - size, cz],
            [cx - size, cy - size, cz],
            [cx - size, cy + size, cz],
            [cx + size, cy + size, cz],
        ],
        dtype=float,
    )


def _plot_xy_trajectory(agents: list[DoubleIntegratorAgent], output_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 6))
    for agent in agents:
        pos_hist = np.array(agent.pos_hist)
        ax.plot(pos_hist[:, 0], pos_hist[:, 1], label=f"agent {agent.id}")
        square = _square_outline(pos_hist[-1, :2], size=UAV_SQUARE_HALF_SIZE)
        ax.plot(square[:, 0], square[:, 1], color="black", linewidth=1.5)
        ax.text(
            pos_hist[-1, 0] + 0.15 * UAV_SQUARE_HALF_SIZE,
            pos_hist[-1, 1] + 0.15 * UAV_SQUARE_HALF_SIZE,
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
        ax.scatter(pos_hist[-1, 0], pos_hist[-1, 1], pos_hist[-1, 2], s=30)
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
    edge_target = np.array([get_desired_distance(t) for t in sim_time[:-1]])

    plt.figure(figsize=(8, 5))
    for name, values in controller.edge_distance_history.items():
        plt.plot(sim_time[:-1], values, label=name)
    plt.plot(sim_time[:-1], edge_target, "--", label="edge target")
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
        plt.plot(sim_time[:-1], np.array(values) - edge_target, label=name)
    plt.axhline(0.0, color="black", linewidth=1.0)
    plt.xlabel("time (s)")
    plt.ylabel("error (m)")
    plt.title("Edge Distance Errors")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "edge_distance_errors.png", dpi=200, bbox_inches="tight")
    plt.close()


def _plot_tensions(
    controller: UndirectedFormationController,
    sim_time: np.ndarray,
    output_dir: Path,
) -> None:
    plt.figure(figsize=(8, 5))
    for name, values in controller.tension_history.items():
        plt.plot(sim_time[:-1], values, label=name)
    plt.xlabel("time (s)")
    plt.ylabel("tension (N)")
    plt.title("Rope Tensions")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "rope_tensions.png", dpi=200, bbox_inches="tight")
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


def _plot_disturbance_estimation(
    agents: list[DoubleIntegratorAgent],
    controller: UndirectedFormationController,
    config: SimulationConfig,
    output_dir: Path,
) -> None:
    if not controller.use_leso:
        return

    time_axis = np.array(controller.control_time_history)
    control_stride = int(round(config.control_dt / config.sim_dt))
    for agent_idx in range(4):
        key = f"agent_{agent_idx}"
        estimates = np.array(controller.dist_hat_history[key])
        if estimates.size == 0:
            continue

        agent = agents[agent_idx]
        rope_force_full = np.array(agent.rope_force_hist)
        step_dist_full = np.array(agent.disturbance_hist)
        actual_external_full = rope_force_full / config.mass + step_dist_full

        sampled_indices = [0]
        if len(actual_external_full) > 1:
            sampled_indices.extend(range(control_stride - 1, len(actual_external_full), control_stride))
        sampled_indices = sampled_indices[: len(time_axis)]
        actual = actual_external_full[sampled_indices]
        estimates = estimates[: len(actual)]
        time_plot = time_axis[: len(actual)]

        fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
        labels = ["x", "y", "z"]
        for axis in range(3):
            axes[axis].plot(time_plot, actual[:, axis], label="actual external")
            axes[axis].plot(time_plot, estimates[:, axis], "--", label="estimated")
            axes[axis].set_ylabel(labels[axis])
            axes[axis].grid(True)
            axes[axis].legend()
        axes[0].set_title(f"Agent {agent_idx + 1} Disturbance Estimation")
        axes[-1].set_xlabel("time (s)")
        plt.tight_layout()
        plt.savefig(output_dir / f"disturbance_agent_{agent_idx + 1}.png", dpi=200, bbox_inches="tight")
        plt.close(fig)


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
    uav_lines = []
    labels = []

    for idx, agent in enumerate(agents):
        color = colors[idx % len(colors)]
        traj_line, = ax.plot([], [], [], linewidth=1.5, color=color, label=f"uav_{agent.id}")
        uav_line, = ax.plot([], [], [], linewidth=2.0, color=color)
        label = ax.text(0.0, 0.0, 0.0, "", color="black", fontsize=10)
        traj_lines.append(traj_line)
        uav_lines.append(uav_line)
        labels.append(label)

    ax.legend()

    def update(frame_idx: int):
        pts = sampled_hist[frame_idx]
        actual_frame = frames[frame_idx]

        for idx in range(len(agents)):
            hist = sampled_hist[: frame_idx + 1, idx, :]
            traj_lines[idx].set_data(hist[:, 0], hist[:, 1])
            traj_lines[idx].set_3d_properties(hist[:, 2])

            square = _square_outline_3d(pts[idx], size=UAV_SQUARE_HALF_SIZE)
            uav_lines[idx].set_data(square[:, 0], square[:, 1])
            uav_lines[idx].set_3d_properties(square[:, 2])

            labels[idx].set_position((pts[idx, 0], pts[idx, 1]))
            labels[idx].set_3d_properties(pts[idx, 2] + 0.5 * UAV_SQUARE_HALF_SIZE)
            labels[idx].set_text(f"uav_{idx + 1}")

        ax.set_title(f"Undirected Normal Control Animation\nt = {actual_frame * sim_dt:.2f} s")
        return traj_lines + uav_lines + labels

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
    _plot_tensions(controller, sim_time, output_dir)
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
    _plot_disturbance_estimation(agents, controller, config, output_dir)

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
