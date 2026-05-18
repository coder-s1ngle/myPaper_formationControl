import numpy as np
from agent.observer import LESO
from utils import utils

class formation_control_V_constraints:
    def __init__(self, use_leso=True):
        self.use_leso = use_leso
        self.K_angle = 3.0#3.5#4.2 #7.0
        self.K_damp = 6.0#5.8#6.2#7.0 #7.4
        self.K_dist = 4.0  #4.6 #2.1 #1.75

        self.K_rope = 74#50
        self.C_rope = 18#18
        self.l0 = 1.0


        self.bp_angle = 0.8#1.1#1.35#1.5
        self.d_p_angle = 0.35
        self.mu_p_angle = 0.3

        self.bp_dist = 1.2#1.5#1.8
        self.d_p_dist = 0.2
        self.mu_p_dist = 0.3#0.25#0.12

        self.bp_damp = 1.5#2.0#2.5
        self.d_d_damp = 0.15
        self.mu_d_damp = 0.15

        self.bp_att = 0.8#1.4
        self.d_p_att = 0.1
        self.mu_p_att = 0.3

        self.bp_cop = 1.0#1.9
        self.d_p_cop = 0.05
        self.mu_p_cop = 0.5

        # #角度项
        # self.bp_angle = 2.0#2.5#3.8#4.938   # 饱和幅值基准
        # self.d_p_angle = 0.5#0.5   # 线性区间阈值（rad）
        # self.mu_p_angle = 0.3  # 可变参数∈[0,1]

        # # 距离控制
        # self.bp_dist = 1.8#1.8#1.5#1.8  # 饱和幅值基准
        # self.d_p_dist = 0.3    # 线性区间阈值（m）
        # self.mu_p_dist = 0.12   # 可变参数∈[0,1]

        # # 阻尼项
        # self.bp_damp = 2.5#1.3#1.157#7.2    # 饱和幅值基准
        # self.d_d_damp = 0.15    # 线性区间阈值（m/s）
        # self.mu_d_damp = 0.15   # 可变参数∈[0,1]

        self.K_att = 2.5#4.5#1.5#0.3 #2.5
        self.undirected_normal_hysteresis = 0.05
        self.normal_branch_sign = 1.0
        
        self.norm_vector_error_hist = []

        self.edge_distance_history = {edge["name"]: [] for edge in [
            {"name": "edge12","i":0, "j":1},
            {"name": "edge23","i":1, "j":2},
            {"name": "edge34","i":2, "j":3},
            {"name": "edge41","i":3, "j":0}
        ]}

        self.Tension_history = {edge["name"]: [] for edge in [
            {"name": "edge12","i":0, "j":1},
            {"name": "edge23","i":1, "j":2},
            {"name": "edge34","i":2, "j":3},
            {"name": "edge41","i":3, "j":0}
        ]}
        
        self.Vol_hist = []

        self.formation_angle_error_history = {angle["name"] : [] for angle in [
            {"name": "e_1", "i": 0, "j": 2, "k": 1, "theta_star": np.pi / 2}, 
            {"name": "e_2",   "i": 1, "j": 0, "k": 2, "theta_star": np.pi / 2},  
            {"name": "e_3", "i": 2, "j": 0, "k": 3, "theta_star": np.pi / 2}, 
            {"name": "e_4",   "i": 3, "j": 0, "k": 2, "theta_star": np.pi / 2} 
        ]}

        self.angle_constraints =  [
            {"name": "e_1", "i": 0, "j": 1, "k": 3, "theta_star": np.pi / 2},
            {"name": "e_2",   "i": 1, "j": 2, "k": 0, "theta_star": np.pi / 2},
            {"name": "e_3", "i": 2, "j": 3, "k": 1, "theta_star": np.pi / 2},
            {"name": "e_4",   "i": 3, "j": 0, "k": 2, "theta_star": np.pi / 2}
        ]
        # 添加增益历史记录
        self.gain_history = {
            "agent_0": {
                "angle_gain": [],
                "distance_gain_12": [],
                "distance_gain_14": [],
                "damp_gain_x": [],
                "damp_gain_y": [],
                "damp_gain_z": []
            },
            "agent_1": {
                "angle_gain": [],
                "damp_gain_x": [],
                "damp_gain_y": [],
                "damp_gain_z": []
            },
            "agent_2": {
                "angle_gain": [],
                "damp_gain_x": [],
                "damp_gain_y": [],
                "damp_gain_z": []
            },
            "agent_3": {
                "angle_gain": [],
                "damp_gain_x": [],
                "damp_gain_y": [],
                "damp_gain_z": []
            }
        }

     
        self.edge_constraints = [
            {"name": "edge12","i":0, "j":1},
            {"name": "edge23","i":1, "j":2},
            {"name": "edge34","i":2, "j":3},
            {"name": "edge41","i":3, "j":0}
        ]

        self.leso_dt = 1/50
        self.leso_omega = 8 #10.0#8.0
        self.pose_angle_limit = np.array([30*np.pi/180, 30*np.pi/180])
        self.leso_agents = []
        self.last_applied_acc_cmds = []
        self.filtered_dist_hat = []

    def _observer_label(self):
        return "LESO"

    def _build_agent_leso(self, agent):
        omega_o = self.leso_omega
        observer_kwargs = {
            "beta1": 3 * omega_o,
            "beta2": 3 * omega_o**2,
            "beta3": omega_o**3,
            "dt": self.leso_dt,
        }
        return {
            "x": LESO(
                initial_p=agent.state.pos[0],
                initial_v=agent.state.vel[0],
                **observer_kwargs,
            ),
            "y": LESO(
                initial_p=agent.state.pos[1],
                initial_v=agent.state.vel[1],
                **observer_kwargs,
            ),
            "z": LESO(
                initial_p=agent.state.pos[2],
                initial_v=agent.state.vel[2],
                **observer_kwargs,
            ),
        }

    def reset_observers_from_agents(self, agent_list):
        if not self.use_leso:
            for agent in agent_list:
                agent.controller.history_dist_hat_acc.clear()
                agent.controller.observer_mode_label = "Observer Off"
            return

        self.leso_agents = []
        self.last_applied_acc_cmds = []
        self.filtered_dist_hat = []
        for agent in agent_list:
            self.leso_agents.append(self._build_agent_leso(agent))
            self.last_applied_acc_cmds.append(None)
            self.filtered_dist_hat.append(np.zeros(3, dtype=float))
            agent.controller.history_dist_hat_acc.clear()
            agent.controller.observer_mode_label = self._observer_label()

    def reset_leso_from_agents(self, agent_list):
        self.reset_observers_from_agents(agent_list)

    def _ensure_leso_initialized(self, agent_list):
        if not self.use_leso:
            return

        if len(self.leso_agents) != len(agent_list):

            self.reset_observers_from_agents(agent_list)

    def _compensate_control_with_leso(self, agent_list, u_raw):
        if not self.use_leso:
            return np.array(u_raw, dtype=float, copy=True)

        self._ensure_leso_initialized(agent_list)
        u_comp = np.array(u_raw, dtype=float, copy=True)

        for idx, agent in enumerate(agent_list):
            leso_axes = self.leso_agents[idx]
            prev_applied_acc_cmd = self.last_applied_acc_cmds[idx]
            if prev_applied_acc_cmd is None:
                # Prime the observer with the first control command instead of assuming zero input.
                prev_applied_acc_cmd = np.array(u_raw[idx], dtype=float, copy=True)

            px, py, pz = agent.state.pos
            _, _, a_dist_hat_x = leso_axes["x"].update(px, prev_applied_acc_cmd[0])
            _, _, a_dist_hat_y = leso_axes["y"].update(py, prev_applied_acc_cmd[1])
            _, _, a_dist_hat_z = leso_axes["z"].update(pz, prev_applied_acc_cmd[2])
            dist_hat_raw = np.array([a_dist_hat_x, a_dist_hat_y, a_dist_hat_z], dtype=float)

            self.filtered_dist_hat[idx] = dist_hat_raw
            agent.controller.history_dist_hat_acc.append(dist_hat_raw.copy())
            u_comp[idx] -= dist_hat_raw
            g = agent.params.g
            psi = utils.quat_to_euler(agent.state.quaternion)[2]
            R_psi = np.array([
                [np.cos(psi), -np.sin(psi)],
                [np.sin(psi),  np.cos(psi)]
            ])
            A_psi = R_psi @ np.array([[0, 1], [-1, 0]])
            Thetah_d = np.linalg.solve(-g * A_psi, u_comp[idx, :2])
            Thetah_d = np.clip(Thetah_d, -self.pose_angle_limit, self.pose_angle_limit)
            applied_acc_xy = -g * A_psi @ Thetah_d

            self.last_applied_acc_cmds[idx] = np.array([
                applied_acc_xy[0],
                applied_acc_xy[1],
                u_comp[idx, 2],
            ])

        return u_comp
        
    def unit(self,vec:np.ndarray,eps:float = 1e-12) -> np.ndarray:
        n = np.linalg.norm(vec)
        if n < eps:
            return np.zeros_like(vec)
        return vec / n


    
    
    def saturation_function_on_K(self,error, error_threshold, bp, mu):
        '''饱和控制 系数设置'''
        if abs(error)>error_threshold:
            return bp * (abs(error) ** (mu - 1))
        else:
            return bp * (abs(error_threshold ** (mu - 1)))
    

    # Old unsigned angle definition kept here for quick fallback.
    def angle_at_i(self,pi: np.ndarray, pj: np.ndarray, pk: np.ndarray, eps: float = 1e-12) -> float:
        """计算点i处，由j-i-k构成的夹角（弧度）"""
        rij = pj - pi
        rik = pk - pi
        nij = np.linalg.norm(rij)
        nik = np.linalg.norm(rik)
        if nij < eps or nik < eps:
            return 0.0

        rij = rij / nij
        rik = rik / nik

        dot = float(np.dot(rij, rik))
        cross_vec = np.cross(rij, rik) 
        cross_norm = np.linalg.norm(cross_vec)
        theta = np.arctan2(cross_norm, dot)               # [0, pi]
        return float(theta)

    def get_desired_distance(self, t: float) -> float:
        """Step changes: 1.0 m for t<60s, then 0.4 m for t≥60s."""
        return 1.0 if t < 60.0 else 0.4

    def get_desired_diag13(self, t: float) -> float:
        return self.get_desired_distance(t) * np.sqrt(2)

    def get_desired_normal(self, v_star, t: float) -> np.ndarray:
        """Step change: [1,0,0] for t<60s, [0,0,-1] for t≥60s."""
        if t >= 60.0:
            return np.array([0.0, 0.0, -1.0])
        return np.array([1.0, 0.0, 0.0])

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


    def control_law_output(self, V_star_stack, agent_list, t, p_target=None):
        u = np.zeros((4,3),dtype=float)
        #angle part
        angle_err = {c["name"]: 0.0 for c in self.angle_constraints}
        angle_gains = {i: 0.0 for i in range(4)}
        P_stack = np.vstack([a.state.pos for a in agent_list])
        V_stack = np.vstack([a.state.vel for a in agent_list])
        for angle in self.angle_constraints:
            i = angle["i"]
            j = angle["j"]
            k = angle["k"]
            theta_star = angle["theta_star"]

            zij = self.unit(P_stack[j] - P_stack[i])
            zik = self.unit(P_stack[k] - P_stack[i])

            theta = self.angle_at_i(P_stack[i], P_stack[j], P_stack[k])
            e_angle = theta - theta_star

            self.formation_angle_error_history[angle["name"]].append(e_angle)
            angle_err[angle["name"]] = e_angle

            K_angle_saturation = self.saturation_function_on_K(e_angle, self.d_p_angle, self.bp_angle, self.mu_p_angle)
            angle_gains[i] = K_angle_saturation

            # Undirected angle control uses the heuristic restoring action only.
            u[i] += -K_angle_saturation * e_angle * (zij + zik)
        
        #distance part
        # distance_13 = np.linalg.norm(P_stack[2] - P_stack[0])
        distance_12 = np.linalg.norm(P_stack[1] - P_stack[0])
        distance_14 = np.linalg.norm(P_stack[3] - P_stack[0])
        z12 = self.unit(P_stack[1] - P_stack[0])
        z14 = self.unit(P_stack[3] - P_stack[0])
        e_distance_12 = distance_12 - self.get_desired_distance(t)
        e_distance_14 = distance_14 - self.get_desired_distance(t)
        K_distance_saturation_12 = self.saturation_function_on_K(e_distance_12,self.d_p_dist,self.bp_dist,self.mu_p_dist)
        K_distance_saturation_14 = self.saturation_function_on_K(e_distance_14,self.d_p_dist,self.bp_dist,self.mu_p_dist)
        u[0] += K_distance_saturation_12 * e_distance_12 * z12
        u[0] += K_distance_saturation_14 * e_distance_14 * z14
        u[1] -= K_distance_saturation_12 * e_distance_12 * z12
        u[3] -= K_distance_saturation_14 * e_distance_14 * z14



        # u[0] += self.K_dist * e_distance_12 * z12
        # u[0] += self.K_dist * e_distance_14 * z14
        # u[1] -= self.K_dist * e_distance_12 * z12
        # u[3] -= self.K_dist * e_distance_14 * z14
        #v_star part
        e_V_stack = V_stack - V_star_stack
        damp_gains = {i: [0.0, 0.0, 0.0] for i in range(4)}
        for i in range(V_stack.shape[0]):
            for j in range(V_stack.shape[1]):
                    e_V = e_V_stack[i][j]
                    K_v_damp_saturation = self.saturation_function_on_K(e_V,self.d_d_damp,self.bp_damp,self.mu_d_damp)
                    damp_gains[i][j] = K_v_damp_saturation
                    u[i][j] += -K_v_damp_saturation * e_V
                    # u[i][j] += -self.K_damp * e_V
        #coplanar
        


        z12 = self.unit(P_stack[1] - P_stack[0])   # 1 -> 2
        z13 = self.unit(P_stack[2] - P_stack[0])   # 1 -> 3
        z14 = self.unit(P_stack[3] - P_stack[0])   # 1 -> 4

        n142 = np.cross(z14, z12)
        # Keep the coplanarity normal consistent with the desired "bag opening" direction.
        n142planar = self.unit(n142)
        # n142planar = n142
        # agent3 相对该平面的“离面误差”
        z31 = self.unit(P_stack[0] - P_stack[2])   # 3 -> 1
        z32 = self.unit(P_stack[1] - P_stack[2])   # 3 -> 2
        z34 = self.unit(P_stack[3] - P_stack[2])   # 3 -> 4
        n234_local = self.unit(np.cross(z32, z34))
        coplanar_error = np.dot(z31, n234_local)
        self.Vol_hist.append(coplanar_error)
        K_cop_sat = self.saturation_function_on_K(coplanar_error, self.d_p_cop, self.bp_cop, self.mu_p_cop)
        u[2] -= K_cop_sat * coplanar_error * n234_local
        # # 把共面修正作用在 agent3 上
        # u[2] -= 8.0 * coplanar_error * n142planar #8.0

        #norm vector part
        
        # z12 = P_stack[1] - P_stack[0]
        # z14 = P_stack[3] - P_stack[0]
        # n_curr = np.cross(z12,z14)
        ndes = self.unit(self.get_desired_normal(V_star_stack[0],t))
        selected_ndes = self.unit(self._select_undirected_normal(n142planar, ndes))
        normal_err = np.cross(n142planar,selected_ndes)
        undirected_normal_error = float(
            np.arccos(np.clip(abs(np.dot(n142planar, ndes)), 0.0, 1.0))
        )
        self.norm_vector_error_hist.append(undirected_normal_error)
        P_stack_new = np.vstack([P_stack[0],P_stack[1],P_stack[3]])
        p_c = np.mean(P_stack, axis=0)
        p_c_new = np.mean(P_stack_new,axis = 0)


        for i in range(4):  
            if i!=2:      
                r_i = P_stack[i] - p_c_new
                # u_att = self.K_att * np.cross(normal_err, r_i)
                # u[i] += u_att
                e_n = np.linalg.norm(normal_err)
                K_att_sat = self.saturation_function_on_K(
                    e_n, self.d_p_att, self.bp_att, self.mu_p_att
                )
                u_att = K_att_sat * np.cross(normal_err, r_i)
                u[i] += u_att

        

        # 记录增益历史
        for i in range(4):
            self.gain_history[f"agent_{i}"]["angle_gain"].append(angle_gains.get(i, 0.0))
            self.gain_history[f"agent_{i}"]["damp_gain_x"].append(damp_gains[i][0])
            self.gain_history[f"agent_{i}"]["damp_gain_y"].append(damp_gains[i][1])
            self.gain_history[f"agent_{i}"]["damp_gain_z"].append(damp_gains[i][2])
            if i == 0:
                self.gain_history[f"agent_{i}"]["distance_gain_12"].append(K_distance_saturation_12)
                self.gain_history[f"agent_{i}"]["distance_gain_14"].append(K_distance_saturation_14)
        
        return self._compensate_control_with_leso(agent_list, u)

    

    def Tension_law_output(self,agent_list):
        Fd = np.zeros((4,3),dtype=float)
        P_stack = np.vstack([a.state.pos for a in agent_list])
        V_stack = np.vstack([a.state.vel for a in agent_list])
        for edge in self.edge_constraints:
            i = edge["i"]
            j = edge["j"]
            T = 0.0
            rij = P_stack[j] - P_stack[i]
            d = np.linalg.norm(rij)
            self.edge_distance_history[edge["name"]].append(d)
            zij = rij / d

            if d>self.l0:
                d_dot = float(np.dot(V_stack[j] - V_stack[i], zij))
                T = self.K_rope * (d - self.l0) + self.C_rope * d_dot
                
                T = max(0.0, T)
                Fd[i] += T*zij
                Fd[j] += -T*zij
            
            self.Tension_history[edge["name"]].append(T)
        return Fd
