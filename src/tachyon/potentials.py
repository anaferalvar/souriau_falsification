"""
potentials.py
=============

Potentiels V(phi) pour le tachyon roulant de Sen (2002).

Convention : V est exprimé en unités de la densité critique aujourd'hui
(rho_crit,0 = 1 dans les unités H0=1, 8piG/3=1 utilisées dans
rolling_tachyon.py). Ainsi V(phi) est directement Omega_DE-like.

Pour chaque potentiel, on fournit :
    - V(phi)         : le potentiel lui-même
    - dlnV_dphi(phi) : la dérivée logarithmique V'/V, qui apparaît
                       directement dans l'équation de mouvement
                       (eq. 4 du brief méthodologique) :
                           phi_ddot/(1-phi_dot^2) + 3H phi_dot
                           + V'(phi)/V(phi) = 0
    - with_V0(V0)    : retourne une copie avec une nouvelle normalisation
                       V0, utilisée par la procédure de "shooting"
                       (recherche de V0 qui donne Omega_DE0 = 1-Omega_m
                       aujourd'hui).
"""

import numpy as np


class ExponentialPotential:
    """
    V(phi) = V0 * exp(-lambda * phi)

    Cas le plus étudié dans la littérature (cf. DESI Dark Secrets,
    arXiv:2502.08876). V'/V = -lambda = constante : la "pente
    logarithmique" de la colline ne dépend pas de la position phi.
    """

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
    """
    V(phi) = V0 * phi^(-n),  phi > 0

    Colline en 1/x^n : la pente logarithmique V'/V = -n/phi diminue
    à mesure que phi croît (la bille ralentit en descendant).
    Comportement typiquement plus proche de Lambda à bas redshift.
    """

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
