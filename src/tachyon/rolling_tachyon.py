import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


def tachyon_rhs(N, y, potential, Omega_m, Omega_r=0.0):
    """
    évalue le système d'équations différentielles pour le champ tachyonique.

    calcule les dérivées de phi et phi_dot par rapport au nombre d'e-folds N = ln(a) en integrant
    l'équation de Friedmann et l'équation de mouvement (Klein-Gordon modifiée).

    ARGS:
        N (float): Nombre d'e-folds actuel (ln(a)).
        y (list): État du système [phi, phi_dot].
        potential (object): Instance du potentiel tachyonique.
        Omega_m (float): Densité de matière actuelle.
        Omega_r (float, optionnel): Densité de radiation actuelle. Par défaut à 0.0.

    RETURNS:
        list: Dérivées [dphi/dN, dphi_dot/dN].
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


def solve_background(potential, Omega_m, phi_init, x_init=0.0, Omega_r=0.0, z_init=50.0, n_points=2000):
    """
    intègre l'évolution cosmique du champ tachyonique depuis un redshift initial.

    utilise la méthode Radau pour résoudre le système de z_init jusqu'à aujourd'hui (N=0). (peut-être il y a mieux?)

    ARGS:
        potential (object): Instance du potentiel avec un V0 fixé.
        Omega_m (float): Densité de matière actuelle.
        phi_init (float): Valeur initiale du champ scalaire.
        x_init (float, optionnel): Vitesse initiale du champ (phi_dot). Par défaut à 0.0 (champ gelé).
        Omega_r (float, optionnel): Densité de radiation. Par défaut à 0.0.
        z_init (float, optionnel): Redshift de départ. Par défaut à 50.0.
        n_points (int, optionnel): Résolution temporelle de la solution. Par défaut à 2000.

    RETURNS:
        OdeResult: Objet SciPy contenant les temps d'évaluation (sol.t) et l'état (sol.y).
    """
    N_init = -np.log(1.0 + z_init)
    N_eval = np.linspace(N_init, 0.0, n_points)

    sol = solve_ivp(
        tachyon_rhs, (N_init, 0.0), [phi_init, x_init],
        args=(potential, Omega_m, Omega_r),
        t_eval=N_eval, method="Radau", rtol=1e-9, atol=1e-11
    )
    return sol


def _omega_de_today(V0, potential_template, Omega_m, phi_init, x_init=0.0, Omega_r=0.0, z_init=50.0):
    """
    calcule la densité d'énergie du tachyon aujourd'hui pour un V0 donné.
    -> fonction utilitaire pour l'algorithme de recherche de racine.
    """
    pot = potential_template.with_V0(V0)
    sol = solve_background(pot, Omega_m, phi_init, x_init, Omega_r=Omega_r, z_init=z_init, n_points=400)
    phi0, x0 = sol.y[:, -1]
    return pot.V(phi0) / np.sqrt(1.0 - x0**2)


def shoot_V0(potential_template, Omega_m, phi_init, x_init=0.0, Omega_r=0.0, z_init=50.0, V0_bracket=(1e-4, 50.0)):
    """
    détermine la normalisation V0 du potentiel pour obtenir un univers plat.

    utilise une recherche de racine (algorithme de Brent) pour trouver la valeur
    de V0 telle que Omega_DE(z=0) corresponde exactement à 1 - Omega_m.

    ARGS:
        potential_template (object): Instance du modèle de potentiel.
        Omega_m (float): Densité de matière cible.
        phi_init (float): Condition initiale du champ.
        x_init (float, optionnel): Vitesse initiale du champ. Par défaut à 0.0.
        Omega_r (float, optionnel): Densité de radiation. Par défaut à 0.0.
        z_init (float, optionnel): Redshift de départ. Par défaut à 50.0.
        V0_bracket (tuple, optionnel): Intervalle de recherche pour V0. Par défaut à (1e-4, 50.0).

    RETURNS:
        float: Valeur de V0 résolvant la contrainte géométrique.
    """
    target = 1.0 - Omega_m

    def f(V0):
        return _omega_de_today(V0, potential_template, Omega_m, phi_init, x_init, Omega_r, z_init) - target

    V0_sol = brentq(f, *V0_bracket, xtol=1e-12, rtol=1e-10)
    return V0_sol


def build_tachyon_background(potential_template, Omega_m, phi_init, x_init=0.0, Omega_r=0.0, z_init=50.0, n_points=2000):
    """
    pipeline complet: calibre le potentiel et génère l'historique cosmologique.

    1. Trouve le bon V0 par méthode de tir.
    2. Intègre la trajectoire complète du champ.
    3. Calcule les grandeurs cosmologiques dérivées (H, w, Omega_DE) sur la grille temporelle.

    RETURNS:
        dict: Dictionnaire contenant les tableaux d'évolution (N, a, z, phi, phi_dot, w, H, Omega_DE) et l'objet potentiel final calibré.
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

def sanity_checks(bg, verbose=True):
    """
    vérifie l'intégrité physique et numérique du fond cosmologique généré.

    S'assure du respect des limites théoriques du tachyon de Sen (causalité,
    équation d'état) et des contraintes d'univers plat aux conditions limites.

    Args:
        bg (dict): Dictionnaire de l'historique cosmologique généré.
        verbose (bool, optionnel): Affiche les résultats dans la console. Par défaut à True.

    Returns:
        dict: Résultats des tests sous forme de booléens.
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
        print("SANITY CHECKS")
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