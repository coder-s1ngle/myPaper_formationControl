"""Generate paper-ready comparison figures for Experiment 1 (coplanarity ablation).

Produces three standalone figures (matching paper Section V-A):
  - fig_exp1_angle.pdf    : four angle errors, (a) with / (b) without coplanarity
  - fig_exp1_coplanar.pdf : coplanarity error comparison (single column)
  - fig_exp1_height.pdf   : agent 3 height comparison (single column)

Initial conditions (paper eq. 43):
  p1(0)=[0,0.6,0]^T, p2(0)=[0,0,0]^T, p3(0)=[0.6,0,0.05]^T, p4(0)=[0.6,0.6,0]^T.
Agent #3 has a deliberate 0.05 m out-of-plane perturbation.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Exp_1.exp_1_with_coplanar import simulate as sim_w, SimulationConfig as CfgW
from Exp_1.exp_1_without_coplanar import simulate as sim_wo, SimulationConfig as CfgWO

# ── IEEE-compatible style ────────────────────────────────────────
IEEE_SINGLE_COL = 3.5   # inches
IEEE_DOUBLE_COL = 7.0

plt.rc("font", size=7, family="serif")
plt.rc("axes", titlesize=8, labelsize=7)
plt.rc("xtick", labelsize=6)
plt.rc("ytick", labelsize=6)
plt.rc("legend", fontsize=6)
plt.rc("lines", linewidth=0.9, markersize=2.5)

COLOR_W  = "#2166AC"   # blue — with coplanar
COLOR_WO = "#B2182B"   # red  — without coplanar
LW = 0.8


def run_sims():
    cfg_w = CfgW(save_animation=False)
    agents_w, ctrl_w, _ = sim_w(cfg_w)

    cfg_wo = CfgWO(save_animation=False)
    agents_wo, ctrl_wo, _ = sim_wo(cfg_wo)

    return ctrl_w, ctrl_wo, agents_w, agents_wo


def _common_ax(ax):
    ax.axhline(0, color="gray", linewidth=0.4, linestyle="--")
    ax.set_xlabel("Time (s)")
    ax.grid(True, linewidth=0.3, alpha=0.6)
    ax.tick_params(direction="in")


def _save(fig, name: str, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ── individual figures ───────────────────────────────────────────

def fig_angle_error(ctrl_w, ctrl_wo, out_dir: Path):
    t_w  = np.array(ctrl_w.control_time_history)
    t_wo = np.array(ctrl_wo.control_time_history)

    names = ["e_1", "e_2", "e_3", "e_4"]
    colors4 = ["#2166AC", "#4DAF4A", "#FF7F00", "#94318E"]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(IEEE_SINGLE_COL, 3.0),
                                   constrained_layout=True)

    for i, name in enumerate(names):
        ax1.plot(t_w, np.array(ctrl_w.angle_error_history[name]),
                 color=colors4[i], linewidth=LW,
                 label=f"$e_{{\\theta_{i+1}}}$")
        ax2.plot(t_wo, np.array(ctrl_wo.angle_error_history[name]),
                 color=colors4[i], linewidth=LW)

    ax1.set_title("(a) With $k_{\\mathrm{cop}}$", loc="left", fontsize=7, fontweight="bold")
    ax2.set_title("(b) Without $k_{\\mathrm{cop}}$", loc="left", fontsize=7, fontweight="bold")
    ax1.set_ylabel("Angle error (rad)")
    for ax in (ax1, ax2):
        ax.axhline(0, color="gray", linewidth=0.4, linestyle="--")
        ax.set_xlabel("Time (s)")
        ax.grid(True, linewidth=0.3, alpha=0.6)
        ax.tick_params(direction="in")
    ax1.legend(frameon=False, ncol=2, fontsize=6, loc="upper right")

    _save(fig, "fig_exp1_angle", out_dir)


def fig_coplanar_error(ctrl_w, ctrl_wo, out_dir: Path):
    t_w  = np.array(ctrl_w.control_time_history)
    t_wo = np.array(ctrl_wo.control_time_history)
    c_w  = np.array(ctrl_w.coplanar_error_history)
    c_wo = np.array(ctrl_wo.coplanar_error_history)

    fig, ax = plt.subplots(figsize=(IEEE_SINGLE_COL, 2.2))
    ax.plot(t_w,  c_w,  color=COLOR_W,  linewidth=LW, label="with $k_{\\mathrm{cop}}$")
    ax.plot(t_wo, c_wo, color=COLOR_WO, linewidth=LW, label="without $k_{\\mathrm{cop}}$")
    _common_ax(ax)
    ax.set_ylabel("$e_{\\mathrm{cop}}$")
    ax.legend(frameon=False, loc="upper right")
    _save(fig, "fig_exp1_coplanar", out_dir)


def fig_height(ctrl_w, ctrl_wo, agents_w, agents_wo, out_dir: Path):
    t_w  = np.array(ctrl_w.control_time_history)
    t_wo = np.array(ctrl_wo.control_time_history)

    # agent 3 = index 2; height = -z (so upward is positive)
    z3_w  = -np.array(agents_w[2].pos_hist)[1:, 2]
    z3_wo = -np.array(agents_wo[2].pos_hist)[1:, 2]
    if len(z3_w)  != len(t_w):   z3_w  = z3_w[: len(t_w)]
    if len(z3_wo) != len(t_wo):  z3_wo = z3_wo[: len(t_wo)]

    fig, ax = plt.subplots(figsize=(IEEE_SINGLE_COL, 2.2))
    ax.plot(t_w,  z3_w,  color=COLOR_W,  linewidth=LW, label="with $k_{\\mathrm{cop}}$")
    ax.plot(t_wo, z3_wo, color=COLOR_WO, linewidth=LW, label="without $k_{\\mathrm{cop}}$")
    _common_ax(ax)
    ax.set_ylabel("Height of agent 3 (m)")
    ax.legend(frameon=False, loc="upper right")
    _save(fig, "fig_exp1_height", out_dir)


# ── main ─────────────────────────────────────────────────────────

if __name__ == "__main__":
    ctrl_w, ctrl_wo, agents_w, agents_wo = run_sims()
    out_dir = ROOT / "output" / "double_integrator" / "exp_1_paper"

    fig_angle_error(ctrl_w, ctrl_wo, out_dir)
    fig_coplanar_error(ctrl_w, ctrl_wo, out_dir)
    fig_height(ctrl_w, ctrl_wo, agents_w, agents_wo, out_dir)

    print(f"3 figures saved to {out_dir}")
