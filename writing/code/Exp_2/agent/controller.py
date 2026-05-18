import numpy as np
from utils import utils


class Controller:
    def __init__(self, agent):
        self.agent = agent

        self.pos_loop_control_dt = 1 / 50
        self.att_loop_control_dt = 1 / 200
        self.inner_loop_times = int(
            self.pos_loop_control_dt / self.att_loop_control_dt
        )

        self.twist_lim_body = [
            np.inf,
            np.inf,
            np.inf,
            200 * np.pi / 180,
            200 * np.pi / 180,
            60 * np.pi / 180,
        ]
        self.accel_lim_body = [
            np.inf,
            np.inf,
            np.inf,
            2000 * np.pi / 180,
            2000 * np.pi / 180,
            600 * np.pi / 180,
        ]
        self.pose_lim_body = [
            np.inf,
            np.inf,
            np.inf,
            30 * np.pi / 180,
            30 * np.pi / 180,
            np.pi,
        ]

        self.k_px = np.array([1.4, 0.0, 0.3])
        self.k_py = np.array([1.4, 0.0, 0.3])
        self.k_pz = np.array([4.5, 0.0, 0.05])

        self.k_vx = np.array([4.0, 0.0, 0.3])
        self.k_vy = np.array([4.0, 0.0, 0.3])
        self.k_vz = np.array([30.0, 0.0, 0.05])

        self.k_phi = np.array([10.0, 0.0, 0.0])
        self.k_theta = np.array([10.0, 0.0, 0.0])
        self.k_psi = np.array([10.0, 0.0, 0.0])

        self.k_p = np.array([30.0, 0.0, 0.01])
        self.k_q = np.array([30.0, 0.0, 0.01])
        self.k_r = np.array([30.0, 0.0, 0.01])

        self.integral_e_vh = np.zeros(2)
        self.last_e_vh = np.zeros(2)

        self.integral_e_vz = 0.0
        self.last_e_vz = 0.0

        self.integral_e_omega = np.zeros(3)
        self.last_e_omega = np.zeros(3)

        # Observer traces are recorded by the formation controller.
        self.history_dist_hat_acc = []
        self.history_f_hat = self.history_dist_hat_acc
        self.observer_mode_label = "Observer"

        # Held references for the multi-rate outer/inner loop update.
        self.last_f_d = self.agent.params.mass_b * self.agent.params.g
        self.last_att_cmd = np.zeros(3)

    def pos_controller(self, control_cmd):
        px, py, pz = self.agent.state.pos
        vx, vy, vz = self.agent.state.vel
        _, _, psi = utils.quat_to_euler(self.agent.state.quaternion)

        m = self.agent.params.mass_b
        g = self.agent.params.g

        px_d = control_cmd[0]
        py_d = control_cmd[1]

        k_ph_p = np.array([self.k_px[0], self.k_py[0]])
        ph = np.array([px, py])
        ph_d = np.array([px_d, py_d])
        vh_d = -k_ph_p * (ph - ph_d)
        vh_d = np.clip(
            vh_d,
            -np.array(self.twist_lim_body[0:2]),
            np.array(self.twist_lim_body[0:2]),
        )

        vh = np.array([vx, vy])
        k_vh_p = np.array([self.k_vx[0], self.k_vy[0]])
        k_vh_i = np.array([self.k_vx[1], self.k_vy[1]])
        k_vh_d = np.array([self.k_vx[2], self.k_vy[2]])

        e_vh = vh - vh_d
        self.integral_e_vh += e_vh * self.pos_loop_control_dt
        integral_e_vh = self.integral_e_vh
        differential_e_vh = (e_vh - self.last_e_vh) / self.pos_loop_control_dt

        r_psi = np.array(
            [
                [np.cos(psi), -np.sin(psi)],
                [np.sin(psi), np.cos(psi)],
            ]
        )
        a_psi = r_psi @ np.array([[0, 1], [-1, 0]])
        d_vh_d = (
            -k_vh_p * e_vh
            - k_vh_i * integral_e_vh
            - k_vh_d * differential_e_vh
        )
        d_vh_d = np.clip(
            d_vh_d,
            -np.array(self.accel_lim_body[0:2]),
            np.array(self.accel_lim_body[0:2]),
        )
        thetah_d = np.linalg.solve(-g * a_psi, d_vh_d)
        thetah_d = np.clip(
            thetah_d,
            -np.array(self.pose_lim_body[3:5]),
            np.array(self.pose_lim_body[3:5]),
        )
        phi_d, theta_d = thetah_d

        pz_d = control_cmd[2]
        k_pz_p = self.k_pz[0]
        k_vz_p = self.k_vz[0]
        k_vz_i = self.k_vz[1]
        k_vz_d = self.k_vz[2]

        vz_d = -k_pz_p * (pz - pz_d)
        vz_d = np.clip(vz_d, -self.twist_lim_body[2], self.twist_lim_body[2])
        e_vz = vz - vz_d
        self.integral_e_vz += e_vz * self.pos_loop_control_dt
        integral_e_vz = self.integral_e_vz
        differential_e_vz = (e_vz - self.last_e_vz) / self.pos_loop_control_dt
        d_vz_d = -k_vz_p * e_vz - k_vz_i * integral_e_vz - k_vz_d * differential_e_vz
        d_vz_d = np.clip(
            d_vz_d,
            -self.accel_lim_body[2],
            self.accel_lim_body[2],
        )

        f_d = m * (g - d_vz_d)

        self.last_e_vh = e_vh
        self.last_e_vz = e_vz

        return f_d, phi_d, theta_d

    def formation_controller_input_u(self, control_cmd, U_d):
        g = self.agent.params.g
        d_vh_d = np.array([U_d[0], U_d[1]])

        psi = utils.quat_to_euler(self.agent.state.quaternion)[2]
        r_psi = np.array(
            [
                [np.cos(psi), -np.sin(psi)],
                [np.sin(psi), np.cos(psi)],
            ]
        )
        a_psi = r_psi @ np.array([[0, 1], [-1, 0]])
        thetah_d = np.linalg.solve(-g * a_psi, d_vh_d)
        thetah_d = np.clip(
            thetah_d,
            -np.array(self.pose_lim_body[3:5]),
            np.array(self.pose_lim_body[3:5]),
        )
        phi_d, theta_d = thetah_d[0], thetah_d[1]

        d_vz_d = U_d[2]
        f_d = self.agent.params.mass_b * (g - d_vz_d)
        return f_d, phi_d, theta_d

    def att_controller(self, control_cmd):
        quat = self.agent.state.quaternion
        phi, theta, psi = utils.quat_to_euler(quat)
        p, q, r = self.agent.state.omega

        k_omega_p = np.array([self.k_p[0], self.k_q[0], self.k_r[0]])
        k_omega_i = np.array([self.k_p[1], self.k_q[1], self.k_r[1]])
        k_omega_d = np.array([self.k_p[2], self.k_q[2], self.k_r[2]])
        k_theta_p = np.array([self.k_phi[0], self.k_theta[0], self.k_psi[0]])

        theta_now = np.array([phi, theta, psi])
        theta_d = np.array([control_cmd[0], control_cmd[1], control_cmd[2]])

        omega_d = -k_theta_p * (theta_now - theta_d)
        w_map = np.array(
            [
                [1, np.sin(phi) * np.tan(theta), np.cos(phi) * np.tan(theta)],
                [0, np.cos(phi), -np.sin(phi)],
                [0, np.sin(phi) / np.cos(theta), np.cos(phi) / np.cos(theta)],
            ]
        )
        omega_d = np.linalg.solve(w_map, omega_d)
        omega_d = np.clip(
            omega_d,
            -np.array(self.twist_lim_body[3:6]),
            np.array(self.twist_lim_body[3:6]),
        )

        omega = np.array([p, q, r])
        e_omega = omega - omega_d
        self.integral_e_omega += e_omega * self.att_loop_control_dt
        integral_e_omega = self.integral_e_omega
        differential_e_omega = (e_omega - self.last_e_omega) / self.att_loop_control_dt
        d_omega_d = (
            -k_omega_p * e_omega
            - k_omega_i * integral_e_omega
            - k_omega_d * differential_e_omega
        )
        d_omega_d = np.clip(
            d_omega_d,
            -np.array(self.accel_lim_body[3:6]),
            np.array(self.accel_lim_body[3:6]),
        )

        tao_d = self.agent.params.Ib @ d_omega_d
        self.last_e_omega = e_omega
        return tao_d[0], tao_d[1], tao_d[2]

    def update_outer_loop_reference(self, control_cmd, formation_flag, U_d):
        """Update held thrust and attitude references without mutating the plant state."""
        if self.agent.cmd_mode == "pos":
            if formation_flag:
                f_d, phi_d, theta_d = self.formation_controller_input_u(control_cmd, U_d)
            else:
                f_d, phi_d, theta_d = self.pos_controller(control_cmd)
            self.last_f_d = f_d
            self.last_att_cmd = np.array([phi_d, theta_d, 0.0], dtype=float)
            return

        if self.agent.cmd_mode == "att":
            m = self.agent.params.mass_b
            g = self.agent.params.g
            self.last_f_d = m * g / (np.cos(control_cmd[0]) * np.cos(control_cmd[1]))
            self.last_att_cmd = np.array(control_cmd, dtype=float)

    def update_attitude_control(self):
        """Run the attitude loop against the held reference and map to motor speeds."""
        taox_d, taoy_d, taoz_d = self.att_controller(self.last_att_cmd)
        return self.allocate_control([self.last_f_d, taox_d, taoy_d, taoz_d])

    def control_cycle(self, control_cmd, formation_flag, U_d):
        """Compatibility wrapper for callers that still update both loops together."""
        self.update_outer_loop_reference(control_cmd, formation_flag, U_d)
        taox_d, taoy_d, taoz_d = self.att_controller(self.last_att_cmd)
        return self.last_f_d, taox_d, taoy_d, taoz_d

    def allocate_control(self, control_cmd):
        f_d, taox_d, taoy_d, taoz_d = control_cmd
        sqrt2_2 = np.sqrt(2) / 2

        c_t = self.agent.params.C_T
        c_m = self.agent.params.C_M
        arm = self.agent.params.arm_length
        alloc = np.array(
            [
                [c_t, c_t, c_t, c_t],
                [-sqrt2_2 * c_t * arm, sqrt2_2 * c_t * arm, sqrt2_2 * c_t * arm, -sqrt2_2 * c_t * arm],
                [sqrt2_2 * c_t * arm, sqrt2_2 * c_t * arm, -sqrt2_2 * c_t * arm, -sqrt2_2 * c_t * arm],
                [-c_m, c_m, -c_m, c_m],
            ]
        )

        w_squared = np.linalg.solve(alloc, np.array([f_d, taox_d, taoy_d, taoz_d]))
        w_squared = np.maximum(w_squared, 0.0)
        w = np.sqrt(w_squared)

        w = np.maximum(w, 0.0)
        return w[0], w[1], w[2], w[3]
