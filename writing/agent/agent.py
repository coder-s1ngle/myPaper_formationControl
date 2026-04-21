import yaml
import numpy as np
from utils import utils  
from .controller import Controller


def rk4_step(ode_func, state, dt, *args):
    """Fixed-step RK4 integrator for the rigid-body simulation."""
    k1 = np.asarray(ode_func(0.0, state, *args), dtype=float)
    k2 = np.asarray(ode_func(0.5 * dt, state + 0.5 * dt * k1, *args), dtype=float)
    k3 = np.asarray(ode_func(0.5 * dt, state + 0.5 * dt * k2, *args), dtype=float)
    k4 = np.asarray(ode_func(dt, state + dt * k3, *args), dtype=float)
    return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

class AgentParameters:
    def __init__(self):
        filepath_uav_args = "uav_params.yaml"
        # 显式指定 encoding='utf-8' 读取文件
        with open(filepath_uav_args, 'r', encoding='utf-8') as file:
            params = yaml.safe_load(file)

        self.num_rotors = params["uav"]["num_rotors"]
        self.g = params["uav"]["g"]                                     
        self.mass_b = params["uav"]["mass_b"]                             

        self.Ib_xx = params["uav"]["Ib_xx"]                              
        self.Ib_yy = params["uav"]["Ib_yy"]                                
        self.Ib_zz = params["uav"]["Ib_zz"]                                 
        self.C_T = params["uav"]["C_T"]                                     
        self.C_M = params["uav"]["C_M"]                                     
        self.z_r = params["uav"]["z_r"]                            
        self.arm_length = params["uav"]["arm_length"]  

        self.Ib = np.array([
            [self.Ib_xx, 0.0, 0.0],
            [0.0, self.Ib_yy, 0.0],
            [0.0, 0.0, self.Ib_zz]
        ])
        

class AgentState:
    def __init__(self):
        self.pos = np.zeros(3)
        self.vel = np.zeros(3)
        self.u = np.zeros(3)
        self.u_formation_control = np.zeros(3)
        self.Fd = np.zeros(3)
        self.vel_dot = np.zeros(3)
        self.omega = np.zeros(3)
        self.omega_dot = np.zeros(3)
        self.quaternion = np.zeros(4)
        self.quaternion[0] =1
        self.quaternion[1] =0
        self.quaternion[2] =0
        self.quaternion[3] =0
        self.motor_rpms = np.zeros(4)
        
class AgentHistoryStamp:
    def __init__(self,t,pos,quat):
        self.t = t
        self.pos = pos.copy()
        self.quat = quat.copy()

class Agent:
    def __init__(self, agent_id=""):
        self.agent_id = agent_id
        self.params = AgentParameters()
        self.state = AgentState()
        self.cmd_mode = None
        
        
        self.controller = Controller(self)
        self.current_time = 0.0
        self.history_stamp = []
        self.history_u = []
        self.history_u_formation_control = []
        self.history_pos = []
        self.history_att = []
        self.history_vel = []
        self.history_acc = []
        self.history_time = []
        self.history_cmd = []
        self.history_Fd = []
        self.history_dist_force_actual = []
        self.history_dist_acc_actual = []
        self.history_f_actual = self.history_dist_acc_actual

    def update_one_sim_step(self, dt,F_d):
        pos = self.state.pos
        vel = self.state.vel
        vel_dot = self.state.vel_dot
        omega = self.state.omega
        omega_dot = self.state.omega_dot
        quaternion = self.state.quaternion
        motor_rpms = self.state.motor_rpms
        initial_state = np.array([
            pos[0],pos[1],pos[2],
            vel[0],vel[1],vel[2],
            quaternion[0],quaternion[1],quaternion[2],quaternion[3],
            omega[0],omega[1],omega[2]
        ])
        self.state.Fd[0] = F_d[0]
        self.state.Fd[1] = F_d[1]
        self.state.Fd[2] = F_d[2]
        dist_force = F_d.copy()
        dist_acc = dist_force / self.params.mass_b
        self.history_dist_force_actual.append(dist_force)
        self.history_dist_acc_actual.append(dist_acc)
        next_state = rk4_step(dynamics_QUAV, initial_state, dt, self, self.params, self.state.Fd)
        x, y, z, vx, vy, vz, q0, q1, q2, q3, p, q, r = next_state

        self.state.pos[0] = x
        self.state.pos[1] = y
        self.state.pos[2] = z



        self.state.vel[0] = vx
        self.state.vel[1] = vy
        self.state.vel[2] = vz

        quat_next = np.array([q0, q1, q2, q3], dtype=float)
        quat_norm = np.linalg.norm(quat_next)
        if quat_norm > 0.0:
            quat_next /= quat_norm
        self.state.quaternion[0] = quat_next[0]
        self.state.quaternion[1] = quat_next[1]
        self.state.quaternion[2] = quat_next[2]
        self.state.quaternion[3] = quat_next[3]

        self.state.omega[0] = p
        self.state.omega[1] = q
        self.state.omega[2] = r

        att_euler = np.zeros(3)
        att_euler = utils.quat_to_euler(self.state.quaternion)
        self.history_vel.append(self.state.vel.copy())
        self.history_pos.append(self.state.pos.copy())
        self.history_att.append(att_euler)
        self.history_time.append(self.current_time)
        self.history_Fd.append(np.linalg.norm(self.state.Fd))
        
        self.history_u.append(self.state.u.copy())
        self.history_u_formation_control.append(self.state.u_formation_control.copy())

    def update_control_input(self,control_cmd = None,formation_flag = False,U_d = np.zeros(3)):
        self.update_outer_loop_reference(control_cmd, formation_flag, U_d)
        self.update_attitude_control()

    def update_outer_loop_reference(self, control_cmd = None, formation_flag = False, U_d = np.zeros(3)):
        self.controller.update_outer_loop_reference(control_cmd, formation_flag, U_d)
        if not formation_flag and control_cmd is not None:
            self.history_cmd.append(control_cmd.copy())

    def update_attitude_control(self):
        Omega_1, Omega_2, Omega_3, Omega_4 = self.controller.update_attitude_control()

        self.state.motor_rpms[0] = Omega_1
        self.state.motor_rpms[1] = Omega_2
        self.state.motor_rpms[2] = Omega_3
        self.state.motor_rpms[3] = Omega_4


    """
    # ================= 四旋翼架构示意 =================
    #   俯视图 (z 向下为正, NED坐标系)
    #
    #                   x+ (前)
    #                       ^
    # 
    #        (2) ↺                     (1) ↻
    #
    #                                               > y+ (右)
    # 
    #        (3) ↻                     (4) ↺
    #
    #
    #
    #   编号说明:  1 = 前右, 2 = 前左, 3 = 后左, 4 = 后右
    #   旋转方向:  1,3 顺时针 (↻) ; 2,4 逆时针 (↺)
    #   坐标系:    x 前, y 右, z 下
    # ================================================
    """
def dynamics_QUAV(t, state, agent, params, Fd):
    
    x, y, z, vx, vy, vz, q0, q1, q2, q3, p, q, r = state
    quat = np.array([q0, q1, q2, q3])

    Reb = utils.quat_to_rotation_matrix(quat)

    phi ,theta ,psi = utils.quat_to_euler(quat)

    mb = agent.params.mass_b
    g = agent.params.g

    L_r = agent.params.arm_length
    z_r = agent.params.z_r
    C_T = agent.params.C_T
    C_M = agent.params.C_M
    Ib = agent.params.Ib
    E3 = np.eye(3)

    rpm1 = agent.state.motor_rpms[0]
    rpm2 = agent.state.motor_rpms[1]
    rpm3 = agent.state.motor_rpms[2]
    rpm4 = agent.state.motor_rpms[3]

    # Use the integrator state so the ODE remains a true function of the state.
    omega = np.array([p, q, r], dtype=float)

    b3 = np.array([0,0,1])
    e3 = np.array([0,0,1])

    fr1 = -(C_T * rpm1**2 * b3)
    fr2 = -(C_T * rpm2**2 * b3)
    fr3 = -(C_T * rpm3**2 * b3)
    fr4 = -(C_T * rpm4**2 * b3)
    # print("fr1,fr2,fr3,fr4",fr1,fr2,fr3,fr4)
    # f = fr1 + fr2 + fr3 + fr4

    Fr1 = Reb @ fr1
    Fr2 = Reb @ fr2
    Fr3 = Reb @ fr3
    Fr4 = Reb @ fr4

    Fr = Fr1 + Fr2 + Fr3 + Fr4
    Fg = mb * g * e3

    rrb1, rrb2, rrb3, rrb4 = np.array([np.sqrt(2)/2*L_r, np.sqrt(2)/2*L_r, z_r]), \
                             np.array([np.sqrt(2)/2*L_r, -np.sqrt(2)/2*L_r, z_r]), \
                             np.array([-np.sqrt(2)/2*L_r, -np.sqrt(2)/2*L_r, z_r]), \
                             np.array([-np.sqrt(2)/2*L_r, np.sqrt(2)/2*L_r, z_r])
    
    Mt = np.cross(rrb1, fr1) + np.cross(rrb2, fr2) + np.cross(rrb3, fr3) + np.cross(rrb4, fr4)
    Mr = (-1)**1*C_M*rpm1**2 * b3 + \
         (-1)**2*C_M*rpm2**2 * b3 + \
         (-1)**3*C_M*rpm3**2 * b3 + \
         (-1)**4*C_M*rpm4**2 * b3
    Md = np.zeros(3)
    # print("F",Fr+Fg+Fd)
    didt_v_ib = (Fr+Fg+Fd)/mb
    '''更新加速度 包含拉力'''
    ax,ay,az = didt_v_ib
    agent.state.u = np.array([ax,ay,az])
    qvx = np.array([
        [0, -q3, q2],
        [q3, 0, -q1],
        [-q2, q1, 0]
    ])
    dq0 = - 0.5 * np.dot(np.array([q1, q2, q3]), omega)
    dq1, dq2, dq3 = 0.5 * (q0 * np.eye(3) + qvx) @ omega

    domega = np.linalg.solve(Ib,(Mt + Mr + Md - np.cross(omega,Ib @ omega)))

    dp, dq, dr = domega

    # agent.state.vel[0] = ax
    # agent.state.vel[1] = ay
    # agent.state.vel[2] = az

    # agent.state.omega_dot[0] = dp
    # agent.state.omega_dot[1] = dq
    # agent.state.omega_dot[2] = dr

    return [vx, vy, vz, ax, ay, az, dq0, dq1, dq2, dq3, dp, dq, dr]


    


