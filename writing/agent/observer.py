class LESO:
    """Full-order third-order linear ESO driven by position measurements."""

    def __init__(self, beta1, beta2, beta3, dt, initial_p=0.0, initial_v=0.0):
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
        e = p_meas - self.x1_hat
        self.x1_hat += self.dt * (self.x2_hat + self.beta1 * e)
        self.x2_hat += self.dt * (self.x3_hat + u_control + self.beta2 * e)
        self.x3_hat += self.dt * (self.beta3 * e)
        return self.x1_hat, self.x2_hat, self.x3_hat

    def get_states(self):
        return self.x1_hat, self.x2_hat, self.x3_hat

    def get_disturbance(self):
        return self.x3_hat
