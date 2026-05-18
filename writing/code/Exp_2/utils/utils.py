from scipy.spatial.transform import Rotation as R
import numpy as np




def quat_to_euler(q_wxyz, degrees: bool = False):
    """
    四元数 [w, x, y, z] -> 欧拉角 (roll, pitch, yaw)
    roll, pitch, yaw 对应 'xyz' 顺序
    """
    q = np.asarray(q_wxyz, dtype=float)
    assert q.shape == (4,)
    # wxyz -> xyzw
    q_xyzw = [q[1], q[2], q[3], q[0]]
    r = R.from_quat(q_xyzw)
    roll, pitch, yaw = r.as_euler('xyz', degrees=degrees)
    return roll, pitch, yaw


def euler_to_quat(roll: float, pitch: float, yaw: float, degrees: bool = False):
    """
    欧拉角 (roll, pitch, yaw) -> 四元数 [w, x, y, z]
    """
    r = R.from_euler('xyz', [roll, pitch, yaw], degrees=degrees)
    q_xyzw = r.as_quat()  # [x, y, z, w]
    q_wxyz = np.array([q_xyzw[3], q_xyzw[0], q_xyzw[1], q_xyzw[2]], dtype=float)
    return q_wxyz

def euler_to_rotation_matrix(phi: float, theta: float, psi: float,
                             body2inertial: bool = True) -> np.ndarray:
    """
    根据欧拉角 (roll=phi, pitch=theta, yaw=psi) 生成旋转矩阵

    参数:
        phi            : 滚转角 (roll, rad)
        theta          : 俯仰角 (pitch, rad)
        psi            : 偏航角 (yaw, rad)
        body2inertial  : True  -> 返回 机体->惯性 的旋转矩阵
                         False -> 返回 惯性->机体 的旋转矩阵 (即转置)

    返回:
        R : (3,3) numpy array
    """

    Rx = np.array([
        [1.0, 0.0, 0.0],
        [0.0, np.cos(phi), -np.sin(phi)],
        [0.0, np.sin(phi), np.cos(phi)]
    ])

    Ry = np.array([
        [np.cos(theta), 0.0, np.sin(theta)],
        [0.0, 1.0, 0.0],
        [-np.sin(theta), 0.0, np.cos(theta)]
    ])

    Rz = np.array([
        [np.cos(psi), -np.sin(psi), 0.0],
        [np.sin(psi),  np.cos(psi), 0.0],
        [0.0,          0.0,         1.0]
    ])

    # 这里保持你原来习惯：Rb2i = Rx * Ry * Rz （对应 'xyz' 顺序）
    Rb2i = Rx @ Ry @ Rz

    if body2inertial:
        return Rb2i
    else:
        return Rb2i.T  # 惯性 -> 机体
    
def quat_to_rotation_matrix(q_wxyz) -> np.ndarray:
    """
    四元数 [w, x, y, z] -> 旋转矩阵 (3x3)
    """
    q = np.asarray(q_wxyz, dtype=float)
    assert q.shape == (4,)
    # SciPy: [x, y, z, w]
    q_xyzw = [q[1], q[2], q[3], q[0]]
    r = R.from_quat(q_xyzw)
    return r.as_matrix()