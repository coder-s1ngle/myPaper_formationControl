"""Generate combined publication-quality figures for Experiment 1 (coplanarity ablation).

Self-contained: runs both simulations and produces:
  - fig_exp1_trajectory.pdf/png: 3-D trajectories, (a) with coplanarity (b) without
  - fig_exp1_errors.pdf/png: error evolution, (a) coplanar err (b) angle errs
                              (c) distance errs (d) normal err norm
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# ── Helpers (duplicated from experiment scripts for self-containment) ──────

def _unit(vec: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    norm = np.linalg.norm(vec)
    if norm < eps:
        return np.zeros_like(vec)
    return vec / norm


def _angle_at_i(pi: np.ndarray, pj: np.ndarray, pk: np.ndarray,
                eps: float = 1e-12) -> float:
    rij = pj - pi;  rik = pk - pi
    nij = np.linalg.norm(rij);  nik = np.linalg.norm(rik)
    if nij < eps or nik < eps:
        return 0.0
    rij = rij / nij;  rik = rik / nik
    dot = float(np.clip(np.dot(rij, rik), -1.0, 1.0))
    cross_norm = float(np.linalg.norm(np.cross(rij, rik)))
    return float(math.atan2(cross_norm, dot))


# ── Paper-matched parameters ───────────────────────────────────────────────

V_STAR = np.array([0.5, 0.0, 0.0])
DESIRED_DIST = 10.0
DESIRED_NORMAL = np.array([0.0, 0.0, 1.0])
TOTAL_TIME = 60.0
SIM_DT = 1.0 / 500.0
CTRL_DT = 1.0 / 50.0

GAINS = {
    "bp_angle": 0.5,   "bp_dist": 0.3,   "bp_damp": 1.0,
    "bp_att": 0.4,     "bp_cop": 0.4,
}
# thresholds / mus (unused when saturation off, but kept for completeness)
GAIN_THRESH = {"angle": 0.35, "dist": 0.2, "damp": 0.15, "att": 0.1, "cop": 0.05}
GAIN_MU = {"angle": 0.3, "dist": 0.3, "damp": 0.15, "att": 0.3, "cop": 0.5}


class DoubleIntegratorAgent:
    def __init__(self, agent_id: int, p0: list[float],
                 v0: list[float] | None = None):
        self.id = agent_id
        self.pos = np.asarray(p0, dtype=float).reshape(3)
        self.vel = (np.zeros(3, dtype=float) if v0 is None
                    else np.asarray(v0, dtype=float).reshape(3))
        self.pos_hist = [self.pos.copy()]
        self.vel_hist = [self.vel.copy()]

    def step(self, dt: float, formation_cmd: np.ndarray) -> None:
        total_acc = np.asarray(formation_cmd, dtype=float).reshape(3)
        self.pos = self.pos + self.vel * dt + 0.5 * total_acc * dt * dt
        self.vel = self.vel + total_acc * dt
        self.pos_hist.append(self.pos.copy())
        self.vel_hist.append(self.vel.copy())


class FormationController:
    """Implements the paper's controller (5-term) with optional coplanarity."""

    def __init__(self, use_coplanar: bool = True):
        self.use_coplanar = use_coplanar
        self.normal_branch_sign = 1.0
        self.hysteresis = 0.05

        self.angle_constraints = [
            {"name": "e_1", "i": 0, "j": 1, "k": 3, "theta_star": math.pi / 2.0},
            {"name": "e_2", "i": 1, "j": 2, "k": 0, "theta_star": math.pi / 2.0},
            {"name": "e_3", "i": 2, "j": 3, "k": 1, "theta_star": math.pi / 2.0},
            {"name": "e_4", "i": 3, "j": 0, "k": 2, "theta_star": math.pi / 2.0},
        ]

        self.ctrl_time: list[float] = []
        self.angle_err: dict[str, list[float]] = {
            c["name"]: [] for c in self.angle_constraints}
        self.edge_dist: dict[str, list[float]] = {
            "edge12": [], "edge23": [], "edge34": [], "edge41": []}
        self.coplanar_err: list[float] = []
        self.normal_err: list[float] = []
        self.normal_align: list[float] = []

    def _select_normal(self, n_cur: np.ndarray, n_des: np.ndarray) -> np.ndarray:
        a = float(np.dot(n_cur, n_des))
        if a > self.hysteresis:
            self.normal_branch_sign = 1.0
        elif a < -self.hysteresis:
            self.normal_branch_sign = -1.0
        return self.normal_branch_sign * n_des

    def control(self, agents: list[DoubleIntegratorAgent], t: float) -> np.ndarray:
        u = np.zeros((4, 3), dtype=float)
        p = np.vstack([a.pos for a in agents])
        v = np.vstack([a.vel for a in agents])
        v_star = np.tile(V_STAR, (4, 1))

        self.ctrl_time.append(t)

        # (1) Angle constraints
        for c in self.angle_constraints:
            zij = _unit(p[c["j"]] - p[c["i"]])
            zik = _unit(p[c["k"]] - p[c["i"]])
            theta = _angle_at_i(p[c["i"]], p[c["j"]], p[c["k"]])
            e = theta - c["theta_star"]
            self.angle_err[c["name"]].append(e)
            u[c["i"]] += -GAINS["bp_angle"] * e * (zij + zik)

        # (2) Distance constraints (edges (1,2) and (1,4))
        d12 = np.linalg.norm(p[1] - p[0])
        d14 = np.linalg.norm(p[3] - p[0])
        d23 = np.linalg.norm(p[2] - p[1])
        d34 = np.linalg.norm(p[3] - p[2])
        self.edge_dist["edge12"].append(d12)
        self.edge_dist["edge23"].append(d23)
        self.edge_dist["edge34"].append(d34)
        self.edge_dist["edge41"].append(d14)

        z12 = _unit(p[1] - p[0])
        z14 = _unit(p[3] - p[0])
        e12 = d12 - DESIRED_DIST
        e14 = d14 - DESIRED_DIST

        u[0] += GAINS["bp_dist"] * e12 * z12
        u[0] += GAINS["bp_dist"] * e14 * z14
        u[1] -= GAINS["bp_dist"] * e12 * z12
        u[3] -= GAINS["bp_dist"] * e14 * z14

        # (3) Velocity damping
        e_v_stack = v - v_star
        for i in range(4):
            u[i] += -GAINS["bp_damp"] * e_v_stack[i]

        # (4) Coplanarity term (agent 3 only)
        z31 = _unit(p[0] - p[2])
        z32 = _unit(p[1] - p[2])
        z34 = _unit(p[3] - p[2])
        n234 = _unit(np.cross(z32, z34))
        e_cop = float(np.dot(z31, n234))
        self.coplanar_err.append(e_cop)

        if self.use_coplanar:
            u[2] -= GAINS["bp_cop"] * e_cop * n234

        # (5) Undirected normal control
        n142 = _unit(np.cross(z14, z12))
        selected = _unit(self._select_normal(n142, DESIRED_NORMAL))
        e_n = np.cross(n142, selected)
        und_ang = math.acos(float(np.clip(abs(np.dot(n142, DESIRED_NORMAL)),
                                          0.0, 1.0)))
        self.normal_err.append(und_ang)
        self.normal_align.append(float(np.dot(n142, selected)))

        p_center = np.mean(np.vstack([p[0], p[1], p[3]]), axis=0)
        for idx in (0, 1, 3):
            r_i = p[idx] - p_center
            u[idx] += GAINS["bp_att"] * np.cross(e_n, r_i)

        return u


def simulate(use_coplanar: bool):
    agents = [
        DoubleIntegratorAgent(1, [0.0, 10.0, 0.0]),
        DoubleIntegratorAgent(2, [0.0, 0.0, 0.0]),
        DoubleIntegratorAgent(3, [10.0, 0.0, 0.5]),
        DoubleIntegratorAgent(4, [10.0, 10.0, 0.0]),
    ]
    ctrl = FormationController(use_coplanar=use_coplanar)

    t = 0.0
    n_steps = int(round(TOTAL_TIME / SIM_DT))
    sim_time = np.zeros(n_steps + 1, dtype=float)
    ctrl_elapsed = 0.0
    u_latest = ctrl.control(agents, t)

    for step in range(n_steps):
        for i, a in enumerate(agents):
            a.step(SIM_DT, u_latest[i])
        t += SIM_DT
        sim_time[step + 1] = t
        ctrl_elapsed += SIM_DT
        if ctrl_elapsed >= CTRL_DT - 1e-12:
            while ctrl_elapsed >= CTRL_DT - 1e-12:
                ctrl_elapsed -= CTRL_DT
            u_latest = ctrl.control(agents, t)

    return agents, ctrl, sim_time


# ── Plotting params ─────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "legend.fontsize": 6.5,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "lines.linewidth": 0.9,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "text.usetex": False,
})

COLOURS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
OUTPUT_DIR = ROOT / "output" / "double_integrator" / "exp_1_paper_figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
PAPER_FIG_DIR = Path(
    r"E:\work\paper_writing\writing\正文编写\myPaper_formationControl\figures")


# ── Run simulations ─────────────────────────────────────────────────────────
print("Running with coplanarity …")
ag_w, ct_w, t_w = simulate(use_coplanar=True)
print("Running without coplanarity …")
ag_wo, ct_wo, t_wo = simulate(use_coplanar=False)

# ── Figure 1: 3-D trajectories side-by-side ─────────────────────────────────
fig = plt.figure(figsize=(7.0, 3.2))

for sub_idx, (agents, label) in enumerate(
    [(ag_w, "With coplanarity term"),
     (ag_wo, "Without coplanarity term ($k_{\\mathrm{cop}}=0$)")],
    start=1,
):
    ax = fig.add_subplot(1, 2, sub_idx, projection="3d")
    for i, agent in enumerate(agents):
        pos = np.array(agent.pos_hist)
        ax.plot(pos[:, 0], pos[:, 1], pos[:, 2],
                color=COLOURS[i], label=f"Agent {agent.id}",
                linewidth=0.7)
        ax.scatter(*pos[0], s=18, color=COLOURS[i], marker="o", zorder=5)
        ax.scatter(*pos[-1], s=28, color=COLOURS[i], marker="s", zorder=5)
    ax.set_xlabel("x (m)", labelpad=3)
    ax.set_ylabel("y (m)", labelpad=3)
    ax.set_zlabel("z (m)", labelpad=3)
    ax.set_title(f"({chr(96 + sub_idx)}) {label}", fontsize=8.5, pad=3)
    ax.legend(loc="upper left", fontsize=5.5)
    ax.view_init(elev=26, azim=-55)
    ax.tick_params(pad=2)

plt.subplots_adjust(wspace=0.28, left=0.02, right=0.98, bottom=0.06, top=0.93)
fig.savefig(OUTPUT_DIR / "fig_exp1_trajectory.png")
fig.savefig(OUTPUT_DIR / "fig_exp1_trajectory.pdf")
plt.close(fig)
print("Saved fig_exp1_trajectory")

# ── Figure 2: Error evolution ───────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.0))
(ax_cop, ax_ang), (ax_dist, ax_norm) = axes

# (a) Coplanarity error
t_cw = np.array(ct_w.ctrl_time)
t_cwo = np.array(ct_wo.ctrl_time)
ax_cop.plot(t_cw, ct_w.coplanar_err, color=COLOURS[0], label="With coplanarity")
ax_cop.plot(t_cwo, ct_wo.coplanar_err, color=COLOURS[1],
            label="Without coplanarity")
ax_cop.axhline(0, color="grey", linewidth=0.5, linestyle="--")
ax_cop.set_xlabel("Time (s)")
ax_cop.set_ylabel("$e_{\\mathrm{cop}}$")
ax_cop.set_title("(a) Coplanarity error")
ax_cop.legend()
ax_cop.grid(True, alpha=0.3)

# (b) Angle errors
for idx, (name, vals) in enumerate(ct_w.angle_err.items()):
    ax_ang.plot(t_cw, vals, linewidth=0.7, color=COLOURS[idx],
                label=f"{name} (w/)")
for idx, (name, vals) in enumerate(ct_wo.angle_err.items()):
    ax_ang.plot(t_cwo, vals, linewidth=0.7, linestyle="--", color=COLOURS[idx],
                label=f"{name} (w/o)")
ax_ang.axhline(0, color="grey", linewidth=0.5, linestyle="--")
ax_ang.set_xlabel("Time (s)")
ax_ang.set_ylabel("$e_{\\theta_i}$ (rad)")
ax_ang.set_title("(b) Angle errors")
handles, labels = ax_ang.get_legend_handles_labels()
ax_ang.legend(handles[:4], labels[:4], ncol=2, fontsize=5.5)
ax_ang.grid(True, alpha=0.3)

# (c) Distance errors
for key, color, lbl in zip(
    ["edge12", "edge41"], COLOURS[:2], ["e_{12}", "e_{14}"]
):
    vals_w = np.array(ct_w.edge_dist[key]) - DESIRED_DIST
    vals_wo = np.array(ct_wo.edge_dist[key]) - DESIRED_DIST
    ax_dist.plot(t_cw, vals_w, linewidth=0.7, color=color,
                 label=f"${lbl}$ (w/)")
    ax_dist.plot(t_cwo, vals_wo, linewidth=0.7, linestyle="--", color=color,
                 label=f"${lbl}$ (w/o)")
ax_dist.axhline(0, color="grey", linewidth=0.5, linestyle="--")
ax_dist.set_xlabel("Time (s)")
ax_dist.set_ylabel("$e_d$ (m)")
ax_dist.set_title("(c) Distance errors")
ax_dist.legend(fontsize=6.5)
ax_dist.grid(True, alpha=0.3)

# (d) Normal error
ax_norm.plot(t_cw, ct_w.normal_err, color=COLOURS[0], label="With coplanarity")
ax_norm.plot(t_cwo, ct_wo.normal_err, color=COLOURS[1],
             label="Without coplanarity")
ax_norm.set_xlabel("Time (s)")
ax_norm.set_ylabel("$\\|e_n\\|$ (rad)")
ax_norm.set_title("(d) Undirected plane error")
ax_norm.legend()
ax_norm.grid(True, alpha=0.3)

plt.tight_layout()
fig.savefig(OUTPUT_DIR / "fig_exp1_errors.png")
fig.savefig(OUTPUT_DIR / "fig_exp1_errors.pdf")
plt.close(fig)
print("Saved fig_exp1_errors")

# ── Copy to paper figures directory ─────────────────────────────────────────
for name in ("fig_exp1_trajectory", "fig_exp1_errors"):
    for ext in (".pdf", ".png"):
        src = OUTPUT_DIR / f"{name}{ext}"
        dst = PAPER_FIG_DIR / f"{name}{ext}"
        if src.exists():
            dst.write_bytes(src.read_bytes())
            print(f"Copied {name}{ext} → figures/")

# ── Print summary stats ─────────────────────────────────────────────────────
def _summary(agents, ctrl, tag):
    p = np.vstack([a.pos for a in agents])
    v = np.vstack([a.vel for a in agents])
    d12 = np.linalg.norm(p[1] - p[0])
    d14 = np.linalg.norm(p[3] - p[0])
    print(f"\n=== {tag} ===")
    print(f"Final positions:\n{p}")
    print(f"Final velocities:\n{v}")
    print(f"d12 err: {d12 - DESIRED_DIST:.5f}  d14 err: {d14 - DESIRED_DIST:.5f}")
    print(f"Angle errs: {{{', '.join(f'{k}: {v[-1]:.6f}' for k,v in ctrl.angle_err.items())}}}")
    print(f"Coplanar err: {ctrl.coplanar_err[-1]:.6f}")
    print(f"Normal err:   {ctrl.normal_err[-1]:.6f}")
    print(f"Agent 3 z:    {p[2,2]:.6f}")

_summary(ag_w, ct_w, "WITH coplanarity")
_summary(ag_wo, ct_wo, "WITHOUT coplanarity")
print("\nDone.")
