from agent.agent import Agent
from agent import controller
from utils import utils
from world.world import World
from formation_control_law.angle_diag12_14_V4_constraints import formation_control_V_constraints
import random
import numpy as np


def _unit(vec: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    norm = np.linalg.norm(vec)
    if norm < eps:
        return np.zeros_like(vec)
    return vec / norm


def _get_capture_point(
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

    axis_u = _unit(axis_u_vec)
    axis_v = _unit(axis_v_vec)
    return (
        net_center
        + offset_ratio_u * half_span_u * axis_u
        + offset_ratio_v * half_span_v * axis_v
    )


class run:
    def __init__(
        self,
        formation_flag_bool: bool,
        use_leso: bool = True,
    ):
        self.total_time = 0.0
        self.total_step = 0
        self.dt = 1/500
        self.formation_flag = formation_flag_bool

        self.formation_control_law = formation_control_V_constraints(
            use_leso=use_leso,
        )
        self.agent_num = 4
        self.pos_loop_control_dt = 1/50
        self.att_loop_control_dt = 1/200
        self.control_freq = int(round(self.pos_loop_control_dt / self.dt))
        self.pos_loop_elapsed = self.pos_loop_control_dt
        self.att_loop_elapsed = self.att_loop_control_dt
        self.latest_u_formation = np.zeros((self.agent_num, 3), dtype=float)

        '''control frequence(position) 50Hz
           control frequence(angle) 200Hz
           tick 500Hz'''
        
        self.world = World()
        
        self.v_star = np.array([0.1, 0.0, 0.0])
        self.V_star_stack = np.tile(self.v_star, (4, 1))
        # Target: starts at [15, 0, 0], head-on toward formation
        self.p_target = np.array([15.0, 0.0, 0.0])
        self.v_target = np.array([-0.2, 0.0, 0.0])
        self.target_history: list = []

        # Capture disturbance config (t==60s target hits the net)
        self.capture_start_time = 60.0
        self.capture_end_time = 65.0
        self.capture_pull_force = 1.0  # N
        self.capture_point_offset = (0.25, -0.15)
        self.mass = 1.5  # kg (must match uav_params.yaml)
        
        for idx in range(self.agent_num):
            agent = Agent()
            # quat = utils.euler_to_quat(0.0,0.0, random.uniform(-1.57, 1.57))
            quat = utils.euler_to_quat(0.0,0.0,0.0)
            agent.state.quaternion = quat
            agent.agent_id = 'uav_' + str(idx+1)
            agent.cmd_mode = "pos"
            # agent.cmd_mode = "att"
            agent.controller.pos_loop_control_dt = self.pos_loop_control_dt
            agent.controller.att_loop_control_dt = self.att_loop_control_dt
            agent.controller.inner_loop_times = int(agent.controller.pos_loop_control_dt / agent.controller.att_loop_control_dt)
            self.world.agent_list.append(agent)
        


        for index,agent in enumerate(self.world.agent_list):
            match index:
                case 0:
                    self.world.agent_list[index].state.pos[0] = 0.0   
                    self.world.agent_list[index].state.pos[1] = 0.4  
                    self.world.agent_list[index].state.pos[2] = 0.0
                
                case 1:
                    self.world.agent_list[index].state.pos[0] = 0.0   
                    self.world.agent_list[index].state.pos[1] = 0.0 
                    self.world.agent_list[index].state.pos[2] = 0.0

                case 2:
                    self.world.agent_list[index].state.pos[0] = 0.4
                    self.world.agent_list[index].state.pos[1] = 0.0
                    self.world.agent_list[index].state.pos[2] = 0.0

                case 3:
                    self.world.agent_list[index].state.pos[0] = 0.4   
                    self.world.agent_list[index].state.pos[1] = 0.4  
                    self.world.agent_list[index].state.pos[2] = 0.0

        if hasattr(self.formation_control_law, "reset_observers_from_agents"):
            self.formation_control_law.reset_observers_from_agents(self.world.agent_list)

    def get_v_star(self, t: float) -> np.ndarray:
        return np.array([0.1, 0.0, 0.0])
    # def get_stack_vec(self):
    #     P = np.vstack([a.state.pos[0:2] for a in self.world.agent_list])
    #     V = np.vstack([a.state.vel[0:2] for a in self.world.agent_list])
    #     return P,V
    def get_disturbance(self, t, agent_index, p_stack):
        if not (self.capture_start_time <= t <= self.capture_end_time):
            return np.array([0.0, 0.0, 0.0])

        capture_point = _get_capture_point(
            p_stack,
            self.capture_point_offset[0],
            self.capture_point_offset[1],
        )
        pull_direction = _unit(capture_point - p_stack[agent_index])
        return (self.capture_pull_force / self.mass) * pull_direction

    def _get_control_cmd(self, agent):
        control_cmd = np.zeros(3, dtype=float)
        if agent.cmd_mode == "pos":
            control_cmd[0] = agent.state.pos[0]
            control_cmd[1] = agent.state.pos[1]
            control_cmd[2] = agent.state.pos[2]
        elif agent.cmd_mode == "att":
            control_cmd[0] = 0.0
            control_cmd[1] = 0.0
            control_cmd[2] = 5 * np.pi / 180
        return control_cmd

    def step_world(self):
        self.total_step += 1
        self.total_time += self.dt
        self.v_star = self.get_v_star(self.total_time)
        self.V_star_stack = np.tile(self.v_star, (4, 1))
        self.p_target += self.v_target * self.dt
        self.target_history.append(self.p_target.copy())
        self.pos_loop_elapsed += self.dt
        self.att_loop_elapsed += self.dt
        F_d_list = self.formation_control_law.Tension_law_output(self.world.agent_list)

        if self.pos_loop_elapsed >= self.pos_loop_control_dt:
            while self.pos_loop_elapsed >= self.pos_loop_control_dt:
                self.pos_loop_elapsed -= self.pos_loop_control_dt

            if self.formation_flag:
                self.latest_u_formation = self.formation_control_law.control_law_output(
                    self.V_star_stack,
                    self.world.agent_list,
                    self.total_time,
                    self.p_target,
                )
            else:
                self.latest_u_formation.fill(0.0)

            for index, agent in enumerate(self.world.agent_list):
                agent.state.u_formation_control = self.latest_u_formation[index]
                control_cmd = self._get_control_cmd(agent)
                agent.update_outer_loop_reference(
                    control_cmd,
                    self.formation_flag,
                    self.latest_u_formation[index],
                )

        if self.att_loop_elapsed >= self.att_loop_control_dt:
            while self.att_loop_elapsed >= self.att_loop_control_dt:
                self.att_loop_elapsed -= self.att_loop_control_dt

            for agent in self.world.agent_list:
                agent.update_attitude_control()

        for index, agent in enumerate(self.world.agent_list):
            agent.current_time = self.total_time
            p_stack = np.vstack([a.state.pos for a in self.world.agent_list])
            agent.update_one_sim_step(
                self.dt,
                F_d_list[index] + self.get_disturbance(self.total_time, index, p_stack),
            )
