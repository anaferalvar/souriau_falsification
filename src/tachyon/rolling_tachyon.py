"""
rolling_tachyon.py
===================

Dynamique de fond du tachyon roulant de Sen (2002), brief méthodologique
v2.1, §1 et eq. (1)-(4).

Unités utilisées (pour simplifier l'intégration numérique) :
    - H0 = 1            (le temps est mesuré en unités de 1/H0)
    - 8 pi G / 3 = 1    (donc rho_crit,0 = H0^2 = 1)
    - V(phi) est donc directement exprimé en unités de rho_crit,0

Variable d'intégration : N = ln(a)  (nombre d'e-folds), avec a=1
aujourd'hui donc N=0 aujourd'hui.

Équations (brief, eq. 4 et texte) :

    Lagrangien :  L = -V(phi) sqrt(1 - phi_dot^2)
    Densité     :  rho_phi = V(phi) / sqrt(1 - phi_dot^2)
    Pression    :  p_phi   = -V(phi) sqrt(1 - phi_dot^2)
    w           :  w = phi_dot^2 - 1   in [-1, 0]

    Friedmann (plat, matiere + radiation + tachyon) :
        H(N)^2 = Omega_r * a^-4 + Omega_m * a^-3 + rho_phi(N)

    Equation de mouvement (forme covariante du brief, eq. 4),
    réécrite en variable N = ln(a), avec d/dt = H d/dN :

        dphi/dN = phi_dot / H
        dphi_dot/dN = -(1 - phi_dot^2) * [ 3*phi_dot + V'(phi)/(V(phi)*H) ]

Le système (phi, phi_dot) est intégré de N_init (haut redshift, champ
"gelé" : phi_dot_init = 0, comportement thawing typique) jusqu'à N=0
(aujourd'hui).

La normalisation V0 du potentiel est déterminée par une procédure de
"shooting" : on cherche V0 tel que Omega_DE(aujourd'hui) = 1 - Omega_m,
ce qui équivaut à H(N=0) = 1 (cohérence avec H0=1).
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


# ────────────────────────────────────────────────────────────────────
# Équation de mouvement
# ────────────────────────────────────────────────────────────────────

def tachyon_rhs(N, y, potential, Omega_m, Omega_r=0.0):
    """
    Second membre du système d'ODE en variable N = ln(a).

    y = [phi, phi_dot]
    Retourne [dphi/dN, dphi_dot/dN]
    """
    phi, x = y

    # Sécurité numérique : phi_dot^2 < 1 toujours (borne physique
    # imposée par sqrt(1-phi_dot^2) dans le lagrangien)
    x = np.clip(x, -0.999999, 0.999999)

    a = np.exp(N)

    Vphi = potential.V(phi)
    rho_tach = Vphi / np.sqrt(1.0 - x**2)

    H2 = Omega_r * a**-4 + Omega_m * a**-3 + rho_tach
    H = np.sqrt(H2)

    dlnV = potential.dlnV_dphi(phi)

    dphi_dN = x / H
    dx_dN = -(1.0 - x**2) * (3.0 * x + dlnV / H)

    return [dphi_dN, dx_dN]


# ────────────────────────────────────────────────────────────────────
# Intégration du fond
# ────────────────────────────────────────────────────────────────────

def solve_background(potential, Omega_m, phi_init, x_init=0.0,
                      Omega_r=0.0, z_init=50.0, n_points=2000):
    """
    Intègre le système (phi, phi_dot) de z_init jusqu'à z=0 (N=0).

    x_init=0 correspond à un champ "gelé" au redshift initial — le
    comportement "thawing" standard pour une quintessence/tachyon
    qui ne commence à rouler qu'aux époques récentes.

    Retourne l'objet solution de scipy (sol.t = N, sol.y = [phi, x]).
    """
    N_init = -np.log(1.0 + z_init)
    N_eval = np.linspace(N_init, 0.0, n_points)

    sol = solve_ivp(
        tachyon_rhs, (N_init, 0.0), [phi_init, x_init],
        args=(potential, Omega_m, Omega_r),
        t_eval=N_eval, method="Radau", rtol=1e-9, atol=1e-11
    )
    return sol


# ────────────────────────────────────────────────────────────────────
# Procédure de shooting sur V0
# ────────────────────────────────────────────────────────────────────

def _omega_de_today(V0, potential_template, Omega_m, phi_init,
                     x_init=0.0, Omega_r=0.0, z_init=50.0):
    """
    Pour une normalisation V0 donnée, intègre le fond et retourne
    rho_tach(aujourd'hui) = V(phi0)/sqrt(1-x0^2).

    Si H(0)=1 est satisfait, ce nombre doit valoir exactement
    Omega_DE0 = 1 - Omega_m (condition de fermeture de Friedmann
    aujourd'hui : Omega_m + Omega_DE0 = 1).
    """
    pot = potential_template.with_V0(V0)
    sol = solve_background(pot, Omega_m, phi_init, x_init,
                            Omega_r=Omega_r, z_init=z_init, n_points=400)
    phi0, x0 = sol.y[:, -1]
    return pot.V(phi0) / np.sqrt(1.0 - x0**2)


def shoot_V0(potential_template, Omega_m, phi_init, x_init=0.0,
              Omega_r=0.0, z_init=50.0, V0_bracket=(1e-4, 50.0)):
    """
    Recherche par dichotomie (brentq) la valeur de V0 telle que
    Omega_DE(aujourd'hui) = 1 - Omega_m.

    Retourne V0_solution.
    """
    target = 1.0 - Omega_m

    def f(V0):
        return _omega_de_today(V0, potential_template, Omega_m,
                                phi_init, x_init, Omega_r, z_init) - target

    V0_sol = brentq(f, *V0_bracket, xtol=1e-12, rtol=1e-10)
    return V0_sol


# ────────────────────────────────────────────────────────────────────
# Interface haut niveau : fond complet + w(a)
# ────────────────────────────────────────────────────────────────────

def build_tachyon_background(potential_template, Omega_m, phi_init,
                               x_init=0.0, Omega_r=0.0, z_init=50.0,
                               n_points=2000):
    """
    Construit le fond cosmologique complet pour un potentiel donné :
      1. trouve V0 par shooting (Omega_DE0 = 1-Omega_m) ;
      2. intègre (phi, phi_dot) de z_init à z=0 avec ce V0 ;
      3. calcule w(N), H(N), Omega_DE(N) sur toute la trajectoire.

    Retourne un dict avec les tableaux N, a, z, phi, phi_dot, w, H,
    Omega_DE, ainsi que le potentiel final (avec V0 calibré).
    """
    V0 = shoot_V0(potential_template, Omega_m, phi_init, x_init,
                  Omega_r, z_init)
    pot = potential_template.with_V0(V0)

    sol = solve_background(pot, Omega_m, phi_init, x_init,
                            Omega_r=Omega_r, z_init=z_init,
                            n_points=n_points)

    N = sol.t
    phi, x = sol.y
    a = np.exp(N)
    z = 1.0 / a - 1.0

    Vphi = pot.V(phi)
    rho_tach = Vphi / np.sqrt(1.0 - x**2)
    H = np.sqrt(Omega_r * a**-4 + Omega_m * a**-3 + rho_tach)
    Omega_DE = rho_tach / H**2

    w = x**2 - 1.0

    return {
        "potential": pot,
        "V0": V0,
        "N": N, "a": a, "z": z,
        "phi": phi, "phi_dot": x,
        "w": w, "H": H, "Omega_DE": Omega_DE,
        "Omega_m": Omega_m,
    }


# ────────────────────────────────────────────────────────────────────
# Vérifications de cohérence ("definition of done")
# ────────────────────────────────────────────────────────────────────

def sanity_checks(bg, verbose=True):
    """
    Vérifie les propriétés physiques attendues d'un fond tachyon :
      - w in [-1, 0] partout
      - H(z=0) = 1 (cohérence avec H0=1)
      - Omega_DE(z=0) ~ 1 - Omega_m
      - phi_dot^2 < 1 partout (sqrt(1-phi_dot^2) bien défini)

    Retourne un dict de booléens.
    """
    w = bg["w"]
    H = bg["H"]
    Omega_DE = bg["Omega_DE"]
    Omega_m = bg["Omega_m"]
    x = bg["phi_dot"]

    checks = {
        "w_in_range": bool(np.all((w >= -1.0 - 1e-6) & (w <= 0.0 + 1e-6))),
        "H_today_is_1": bool(np.isclose(H[-1], 1.0, atol=1e-4)),
        "Omega_DE_today_correct": bool(
            np.isclose(Omega_DE[-1], 1.0 - Omega_m, atol=1e-3)
        ),
        "phidot_sub_1": bool(np.all(np.abs(x) < 1.0)),
    }

    if verbose:
        print("── Sanity checks ─────────────────────────────")
        for k, v in checks.items():
            status = "OK" if v else "ÉCHEC"
            print(f"  {k:<28} : {status}")
        print(f"  w(z=0)                       = {w[-1]:.4f}")
        print(f"  w range observé              = "
              f"[{w.min():.4f}, {w.max():.4f}]")
        print(f"  H(z=0)                       = {H[-1]:.6f}")
        print(f"  Omega_DE(z=0)                = {Omega_DE[-1]:.4f}  "
              f"(attendu {1 - Omega_m:.4f})")
        print(f"  max|phi_dot|                 = {np.max(np.abs(x)):.4f}")

    return checks
