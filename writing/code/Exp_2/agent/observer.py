import numpy as np


def _leso_derivatives(x1, x2, x3, u, y, b1, b2, b3):
    """Continuous-time LESO state derivatives (per-axis)."""
    e = y - x1
    dx1 = x2 + b1 * e
    dx2 = x3 + u + b2 * e
    dx3 = b3 * e
    return dx1, dx2, dx3


class LESO:
    """Third-order linear ESO integrated with RK4 (per-axis).

    Continuous-time dynamics (bandwidth parameterisation):
        ż₁ = z₂ + β₁ (y − z₁)
        ż₂ = z₃ + u + β₂ (y − z₁)
        ż₃ =      β₃ (y − z₁)
    where  β₁ = 3 ωₒ,  β₂ = 3 ωₒ²,  β₃ = ωₒ³.

    The RK4 integrator treats the measurement *y* and control input *u* as
    zero-order-held over the time step *dt*, giving a 4th-order accurate
    discrete-time approximation of the true solution.
    """

    def __init__(self, beta1, beta2, beta3, dt,
                 initial_p=0.0, initial_v=0.0):
        self.beta1 = beta1
        self.beta2 = beta2
        self.beta3 = beta3
        self.dt = dt

        self.x1_hat = initial_p
        self.x2_hat = initial_v
        self.x3_hat = 0.0

    def reset(self, initial_p=0.0, initial_v=0.0, initial_f=0.0):
        self.x1_hat = initial_p
        self.x2_hat = initial_v
        self.x3_hat = initial_f

    def update(self, p_meas, u_control):
        """One RK4 step with ZOH on p_meas and u_control."""
        x1, x2, x3 = self.x1_hat, self.x2_hat, self.x3_hat
        b1, b2, b3 = self.beta1, self.beta2, self.beta3
        dt = self.dt
        y = p_meas
        u = u_control

        # k1
        k1_1, k1_2, k1_3 = _leso_derivatives(x1, x2, x3, u, y, b1, b2, b3)

        # k2
        k2_1, k2_2, k2_3 = _leso_derivatives(
            x1 + 0.5 * dt * k1_1,
            x2 + 0.5 * dt * k1_2,
            x3 + 0.5 * dt * k1_3,
            u, y, b1, b2, b3)

        # k3
        k3_1, k3_2, k3_3 = _leso_derivatives(
            x1 + 0.5 * dt * k2_1,
            x2 + 0.5 * dt * k2_2,
            x3 + 0.5 * dt * k2_3,
            u, y, b1, b2, b3)

        # k4
        k4_1, k4_2, k4_3 = _leso_derivatives(
            x1 + dt * k3_1,
            x2 + dt * k3_2,
            x3 + dt * k3_3,
            u, y, b1, b2, b3)

        self.x1_hat += (dt / 6.0) * (k1_1 + 2.0 * k2_1 + 2.0 * k3_1 + k4_1)
        self.x2_hat += (dt / 6.0) * (k1_2 + 2.0 * k2_2 + 2.0 * k3_2 + k4_2)
        self.x3_hat += (dt / 6.0) * (k1_3 + 2.0 * k2_3 + 2.0 * k3_3 + k4_3)

        return self.x1_hat, self.x2_hat, self.x3_hat

    def get_states(self):
        return self.x1_hat, self.x2_hat, self.x3_hat

    def get_disturbance(self):
        return self.x3_hat
