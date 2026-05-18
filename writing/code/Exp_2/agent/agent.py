from pathlib import Path

import yaml
import numpy as np
from numba import njit
from utils import utils
from .controller import Controller


# ═══════════════════════════════════════════════
#  Numba-compiled QUAV dynamics & integrator
# ═══════════════════════════════════════════════

@njit
def _quat_to_rot(q0, q1, q2, q3):
    """Numba-compatible quaternion → rotation matrix (3×3)."""
    R = np.empty((3, 3))
    R[0, 0] = q0*q0 + q1*q1 - q2*q2 - q3*q3
    R[0, 1] = 2.0 * (q1*q2 - q0*q3)
    R[0, 2] = 2.0 * (q0*q2 + q1*q3)
    R[1, 0] = 2.0 * (q1*q2 + q0*q3)
    R[1, 1] = q0*q0 - q1*q1 + q2*q2 - q3*q3
    R[1, 2] = 2.0 * (q2*q3 - q0*q1)
    R[2, 0] = 2.0 * (q1*q3 - q0*q2)
    R[2, 1] = 2.0 * (q0*q1 + q2*q3)
    R[2, 2] = q0*q0 - q1*q1 - q2*q2 + q3*q3
    return R


@njit
def _dynamics_quav(state, motor_rpms, mb, g, L_r, z_r, C_T, C_M,
                   Ib_xx, Ib_yy, Ib_zz, Fd):
    """Numba-compiled QUAV rigid-body dynamics.  Returns 13-element derivative."""
    x, y, z, vx, vy, vz, q0, q1, q2, q3, p_val, q_val, r_val = state

    Reb = _quat_to_rot(q0, q1, q2, q3)

    b3 = np.array([0.0, 0.0, 1.0])
    e3 = np.array([0.0, 0.0, 1.0])

    rpm1, rpm2, rpm3, rpm4 = motor_rpms[0], motor_rpms[1], motor_rpms[2], motor_rpms[3]

    fr1 = -(C_T * rpm1 * rpm1) * b3
    fr2 = -(C_T * rpm2 * rpm2) * b3
    fr3 = -(C_T * rpm3 * rpm3) * b3
    fr4 = -(C_T * rpm4 * rpm4) * b3

    Fr1 = Reb @ fr1
    Fr2 = Reb @ fr2
    Fr3 = Reb @ fr3
    Fr4 = Reb @ fr4

    Fr = Fr1 + Fr2 + Fr3 + Fr4
    Fg = mb * g * e3

    s2 = 0.7071067811865476  # sqrt(2)/2
    rrb1 = np.array([ s2*L_r,  s2*L_r, z_r])
    rrb2 = np.array([ s2*L_r, -s2*L_r, z_r])
    rrb3 = np.array([-s2*L_r, -s2*L_r, z_r])
    rrb4 = np.array([-s2*L_r,  s2*L_r, z_r])

    Mt = (np.cross(rrb1, fr1) + np.cross(rrb2, fr2) +
          np.cross(rrb3, fr3) + np.cross(rrb4, fr4))

    Mr = ((-1.0) * C_M * rpm1 * rpm1 * b3 +
          C_M * rpm2 * rpm2 * b3 +
          (-1.0) * C_M * rpm3 * rpm3 * b3 +
          C_M * rpm4 * rpm4 * b3)

    acc = (Fr + Fg + Fd) / mb
    ax, ay, az = acc[0], acc[1], acc[2]

    # quaternion kinematics (q_dot = 0.5 * omega_quat ⊗ q)
    dq0 = -0.5 * (q1 * p_val + q2 * q_val + q3 * r_val)
    dq1 =  0.5 * (q0 * p_val + q2 * r_val - q3 * q_val)
    dq2 =  0.5 * (q0 * q_val + q3 * p_val - q1 * r_val)
    dq3 =  0.5 * (q0 * r_val + q1 * q_val - q2 * p_val)

    # angular acceleration  Ib·ω̇ = Mt + Mr − ω×(Ib·ω)
    Ix, Iy, Iz = Ib_xx, Ib_yy, Ib_zz
    omega = np.array([p_val, q_val, r_val])
    I_omega = np.array([Ix * p_val, Iy * q_val, Iz * r_val])
    cross_term = np.cross(omega, I_omega)
    rhs = Mt + Mr - cross_term
    dp = rhs[0] / Ix
    dq = rhs[1] / Iy
    dr = rhs[2] / Iz

    return np.array([vx, vy, vz, ax, ay, az,
                     dq0, dq1, dq2, dq3,
                     dp, dq, dr])


@njit
def _rk4_step(state, dt, motor_rpms, mb, g, L_r, z_r, C_T, C_M,
              Ib_xx, Ib_yy, Ib_zz, Fd):
    """Fixed-step RK4 integrator (numba-compiled)."""
    k1 = _dynamics_quav(state, motor_rpms, mb, g, L_r, z_r, C_T, C_M,
                        Ib_xx, Ib_yy, Ib_zz, Fd)
    k2 = _dynamics_quav(state + 0.5*dt*k1, motor_rpms, mb, g, L_r, z_r, C_T, C_M,
                        Ib_xx, Ib_yy, Ib_zz, Fd)
    k3 = _dynamics_quav(state + 0.5*dt*k2, motor_rpms, mb, g, L_r, z_r, C_T, C_M,
                        Ib_xx, Ib_yy, Ib_zz, Fd)
    k4 = _dynamics_quav(state + dt*k3, motor_rpms, mb, g, L_r, z_r, C_T, C_M,
                        Ib_xx, Ib_yy, Ib_zz, Fd)
    return state + (dt / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4)


# ═══════════════════════════════════════════════
#  Agent classes (unchanged logic, precomputes Ib_inv)
# ═══════════════════════════════════════════════

class AgentParameters:
    def __init__(self):
        filepath_uav_args = str(Path(__file__).parent.parent / "uav_params.yaml")
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
        self.quaternion[0] = 1
        self.quaternion[1] = 0
        self.quaternion[2] = 0
        self.quaternion[3] = 0
        self.motor_rpms = np.zeros(4)


class AgentHistoryStamp:
    def __init__(self, t, pos, quat):
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

    def update_one_sim_step(self, dt, F_d):
        pos = self.state.pos
        vel = self.state.vel
        quaternion = self.state.quaternion
        omega = self.state.omega

        initial_state = np.array([
            pos[0], pos[1], pos[2],
            vel[0], vel[1], vel[2],
            quaternion[0], quaternion[1], quaternion[2], quaternion[3],
            omega[0], omega[1], omega[2]
        ])

        self.state.Fd[0] = F_d[0]
        self.state.Fd[1] = F_d[1]
        self.state.Fd[2] = F_d[2]
        dist_force = F_d.copy()
        dist_acc = dist_force / self.params.mass_b
        self.history_dist_force_actual.append(dist_force)
        self.history_dist_acc_actual.append(dist_acc)

        p = self.params
        next_state = _rk4_step(
            initial_state, dt,
            self.state.motor_rpms,
            p.mass_b, p.g, p.arm_length, p.z_r, p.C_T, p.C_M,
            p.Ib_xx, p.Ib_yy, p.Ib_zz,
            F_d,
        )

        x, y, z, vx, vy, vz, q0, q1, q2, q3, p_val, q_val, r_val = next_state

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

        self.state.omega[0] = p_val
        self.state.omega[1] = q_val
        self.state.omega[2] = r_val

        att_euler = utils.quat_to_euler(self.state.quaternion)
        self.history_vel.append(self.state.vel.copy())
        self.history_pos.append(self.state.pos.copy())
        self.history_att.append(att_euler)
        self.history_time.append(self.current_time)
        self.history_Fd.append(np.linalg.norm(self.state.Fd))

        # record total acceleration (including gravity + thrust)
        # re-derive from the dynamics for logging
        state_deriv = _dynamics_quav(
            next_state,
            self.state.motor_rpms,
            p.mass_b, p.g, p.arm_length, p.z_r, p.C_T, p.C_M,
            p.Ib_xx, p.Ib_yy, p.Ib_zz,
            F_d,
        )
        self.state.u = np.array([state_deriv[3], state_deriv[4], state_deriv[5]])
        self.history_u.append(self.state.u.copy())
        self.history_u_formation_control.append(self.state.u_formation_control.copy())

    def update_control_input(self, control_cmd=None, formation_flag=False, U_d=np.zeros(3)):
        self.update_outer_loop_reference(control_cmd, formation_flag, U_d)
        self.update_attitude_control()

    def update_outer_loop_reference(self, control_cmd=None, formation_flag=False, U_d=np.zeros(3)):
        self.controller.update_outer_loop_reference(control_cmd, formation_flag, U_d)
        if not formation_flag and control_cmd is not None:
            self.history_cmd.append(control_cmd.copy())

    def update_attitude_control(self):
        Omega_1, Omega_2, Omega_3, Omega_4 = self.controller.update_attitude_control()
        self.state.motor_rpms[0] = Omega_1
        self.state.motor_rpms[1] = Omega_2
        self.state.motor_rpms[2] = Omega_3
        self.state.motor_rpms[3] = Omega_4
