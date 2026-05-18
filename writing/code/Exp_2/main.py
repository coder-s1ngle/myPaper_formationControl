import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from run import run
from utils import utils
from trajectory_plot.trajectory_plot import *


def main():
    formation_flag = True
    use_leso = True
    set_output_dir(r"output/exp_2_target_capture")

    runner = run.run(
        formation_flag_bool=formation_flag,
        use_leso=use_leso,
    )

    while runner.total_time < 90.0:
        runner.step_world()

    plot_disturbance(runner.world.agent_list, runner.dt)
    plot_disturbance_by_axis(runner.world.agent_list, runner.dt)
    plot_xy_trajectory(runner.world.agent_list)
    plot_3d_trajectory(runner.world.agent_list,
                       t_snapshots=[0, 10, 20, 40, 60, 70, 90])
    plot_z_over_time(runner.world.agent_list)
    plot_attitude(runner.world.agent_list)

    plot_formation_edge_err(runner.formation_control_law.edge_distance_history, runner.dt)
    plot_formation_angle_err(runner.formation_control_law.formation_angle_error_history, runner.dt * runner.control_freq)
    plot_formation_tension(runner.formation_control_law.Tension_history, runner.dt)
    plot_Fd(runner.world.agent_list)
    plot_u(runner.world.agent_list)
    plot_u_formation_control(runner.world.agent_list)
    plot_gain_evolution(runner.formation_control_law.gain_history, runner.dt * runner.control_freq)
    plot_vol_error_evolution(runner.formation_control_law.Vol_hist, runner.dt * runner.control_freq)
    plot_norm_error_evolution(runner.formation_control_law.norm_vector_error_hist, runner.dt * runner.control_freq)


if __name__ == "__main__":
    main()
