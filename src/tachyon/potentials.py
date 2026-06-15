import numpy as np


class ExponentialPotential:
    name = "exponential"

    def __init__(self, V0, lam):
        self.V0 = V0
        self.lam = lam

    def V(self, phi):
        return self.V0 * np.exp(-self.lam * phi)

    def dlnV_dphi(self, phi):
        phi = np.asarray(phi, dtype=float)
        return -self.lam * np.ones_like(phi)

    def with_V0(self, V0):
        return ExponentialPotential(V0, self.lam)

    def __repr__(self):
        return f"ExponentialPotential(V0={self.V0:.4g}, lambda={self.lam:.4g})"


class InversePowerPotential:
    name = "inverse_power"

    def __init__(self, V0, n):
        self.V0 = V0
        self.n = n

    def V(self, phi):
        return self.V0 * np.power(phi, -self.n)

    def dlnV_dphi(self, phi):
        return -self.n / phi

    def with_V0(self, V0):
        return InversePowerPotential(V0, self.n)

    def __repr__(self):
        return f"InversePowerPotential(V0={self.V0:.4g}, n={self.n:.4g})"
