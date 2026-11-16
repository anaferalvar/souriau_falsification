"""
cobaya_interface.py (v4 — correctif d'injection)
=================================================

POURQUOI CETTE VERSION
----------------------
Les versions v1-v3 surchargeaient une méthode ``set_camb_params`` qui
n'existe pas dans le wrapper CAMB de cobaya (vérifié sur cobaya 3.6.2) :
la surcharge n'était jamais appelée, w(a) n'était jamais injecté, et
toutes les chaînes « tachyon » et « reconstruction libre » tournaient en
réalité en ΛCDM (paramètres DE = prior pur, chi2 indépendant de lambda).

COMMENT C'EST CORRIGÉ
---------------------
Dans cobaya, le fond cosmologique et les fonctions de transfert sont
calculés par un *helper theory*, ``CambTransfers``, qui appelle
``CAMB.set(params, state)`` puis CAMB. Les paramètres de fond y sont
routés automatiquement. On :

  1. sous-classe ``CambTransfers`` pour qu'il *déclare supporter* les
     paramètres DE (tachyon_lambda, w_bin1, ...) : cobaya les lui route ;
  2. dans son ``calculate``, on retire ces paramètres, on laisse cobaya
     construire le ``CAMBparams`` standard, puis on remplace
     ``CAMBparams.DarkEnergy`` par la table w(a) du modèle ;
  3. on sous-classe ``CAMB`` uniquement pour brancher ce helper.

Test de non-régression : ``tests/test_injection.py`` vérifie que le chi2
BAO+SNe change quand lambda (ou un bin w_i) change.

MODÈLES (option ``de_model`` dans le YAML, ou classes-alias ci-dessous)
---------------------------------------------------------------------
  exponential    : V = V0 exp(-lambda phi)   params: tachyon_lambda, tachyon_phi_init
  inverse_power  : V = V0 phi^-n            params: tachyon_n,      tachyon_phi_init
  freeform       : w(z) en 4 bins           params: w_bin1..w_bin4

Choix de l'objet DE dans CAMB :
  - tachyon  : DarkEnergyFluid (cs2 = 1), w ∈ [-1, 0) par construction ;
  - freeform : DarkEnergyPPF, qui autorise le franchissement de w = -1.
    (Avec DarkEnergyFluid, CAMB rejette les w(a) qui traversent -1 : avec
    stop_at_error: False ces points auraient eu une vraisemblance nulle,
    ce qui aurait *biaisé le test contre le fantôme*.)
  - Les marches du freeform sont adoucies par une tanh de largeur
    ``step_width`` en z (0.02 par défaut) et échantillonnées finement :
    CAMB interpole par spline cubique, une marche brute provoque des
    oscillations parasites (dépassements sous -1 inventés).
"""

import numpy as np
from camb import dark_energy as camb_de
from cobaya.theories.camb import CAMB
from cobaya.theories.camb.camb import CambTransfers

from src.tachyon.potentials import ExponentialPotential, InversePowerPotential
from src.tachyon.rolling_tachyon import build_tachyon_background

Z_BINS = (0.0, 0.3, 0.6, 1.0, 1.5)     # figé — amendement-01 §3

DE_PARAMS = {
    "exponential":   ("tachyon_lambda", "tachyon_phi_init"),
    "inverse_power": ("tachyon_n", "tachyon_phi_init"),
    "freeform":      ("w_bin1", "w_bin2", "w_bin3", "w_bin4"),
    "lcdm":          (),
}
DE_DERIVED = ("w_today", "phidot2_max")


# ─────────────────────────────────────────────────────────────────────
# Tables w(a)
# ─────────────────────────────────────────────────────────────────────

def freeform_w_table(w_bins, z_bins=Z_BINS, step_width=0.02, z_max=1100.0,
                     n_points=600):
    """w(a) en escalier adouci. w = -1 au-delà de z_bins[-1]."""
    z = np.concatenate([np.linspace(0.0, 2.5, n_points),
                        np.geomspace(2.6, z_max, 60)])
    w_levels = list(w_bins) + [-1.0]
    w = np.full_like(z, w_levels[0])
    for k, edge in enumerate(z_bins[1:]):          # transitions aux bornes
        jump = w_levels[k + 1] - w_levels[k]
        w = w + jump * 0.5 * (1.0 + np.tanh((z - edge) / step_width))
    a = 1.0 / (1.0 + z)
    return a[::-1], w[::-1]              # CAMB exige a croissant, finissant à a = 1


def tachyon_w_table(potential, Omega_m, phi_init, Omega_r, z_init, n_points):
    """Résout le fond tachyon et renvoie (a, w, w_today, phidot2_max)."""
    bg = build_tachyon_background(potential, Omega_m, phi_init=phi_init,
                                  x_init=0.0, Omega_r=Omega_r,
                                  z_init=z_init, n_points=n_points)
    a, w = bg["a"], bg["w"]
    if not np.all(np.isfinite(w)) or np.any(w < -1 - 1e-6) \
            or np.any(np.abs(bg["phi_dot"]) >= 1.0):
        return None
    # au-delà de z_init : champ gelé, w = -1
    a_early = np.geomspace(1e-4, a[0] * 0.999, 40)
    a = np.concatenate([a_early, a])     # CAMB exige une table finissant à a = 1
    w = np.concatenate([np.full(40, -1.0), w])
    return a, w, float(w[-1]), float(np.max(bg["phi_dot"] ** 2))


# ─────────────────────────────────────────────────────────────────────
# Helper theory : c'est ici que w(a) est réellement injecté
# ─────────────────────────────────────────────────────────────────────

class DETransfers(CambTransfers):

    def get_can_support_params(self):
        return set(super().get_can_support_params()) | set(
            DE_PARAMS[self.cobaya_camb.de_model])

    def get_can_provide_params(self):
        return list(DE_DERIVED)

    def _dark_energy(self, de_vals, camb_params):
        """Construit l'objet DarkEnergy CAMB. None => vraisemblance nulle."""
        model = self.cobaya_camb.de_model
        h = camb_params.H0 / 100.0
        Omega_m = (camb_params.omch2 + camb_params.ombh2 + camb_params.omnuh2) / h**2

        if model == "lcdm":
            return None, -1.0, 0.0

        if model == "freeform":
            w_bins = [de_vals[p] for p in DE_PARAMS["freeform"]]
            a, w = freeform_w_table(w_bins, step_width=self.cobaya_camb.step_width)
            de = camb_de.DarkEnergyPPF()
            de.set_w_a_table(a, w)
            return de, float(w_bins[0]), float("nan")

        phi_init = de_vals["tachyon_phi_init"]
        if model == "exponential":
            pot = ExponentialPotential(V0=1.0, lam=de_vals["tachyon_lambda"])
        else:
            pot = InversePowerPotential(V0=1.0, n=de_vals["tachyon_n"])
        c = self.cobaya_camb
        out = tachyon_w_table(pot, Omega_m, phi_init, c.Omega_r, c.z_init, c.n_points)
        if out is None:
            return False, None, None
        a, w, w_today, phidot2_max = out
        de = camb_de.DarkEnergyFluid(cs2=1.0)
        de.set_w_a_table(a, w)
        return de, w_today, phidot2_max

    def calculate(self, state, want_derived=True, **params_values_dict):
        names = DE_PARAMS[self.cobaya_camb.de_model]
        de_vals = {p: params_values_dict.pop(p) for p in names}

        camb_params = self.cobaya_camb.set(params_values_dict, state)
        if not camb_params:
            return False
        try:
            de, w_today, phidot2_max = self._dark_energy(de_vals, camb_params)
        except Exception as e:                       # échec du tir, EDO...
            if self.stop_at_error:
                raise
            self.log.debug("Fond DE échoué : %r", e)
            return False
        if de is False:
            return False
        if de is not None:
            camb_params.DarkEnergy = de
        state["derived"] = {"w_today": w_today, "phidot2_max": phidot2_max}

        try:
            if self.non_linear_sources:
                results = self.camb.get_transfer_functions(
                    camb_params, only_time_sources=True)
            else:
                results = (self.camb.get_transfer_functions(camb_params)
                           if self.needs_perts else self.camb.get_background(camb_params))
            state["results"] = (camb_params, results)
        except self.camb.baseconfig.CAMBError as e:
            if self.stop_at_error:
                raise
            self.log.debug("Calcul CAMB échoué : %s", e)
            return False


class DarkEnergyCAMB(CAMB):
    """CAMB avec énergie noire tabulée. Option YAML : de_model."""

    de_model: str = "exponential"
    Omega_r: float = 9.15e-5
    z_init: float = 50.0
    n_points: int = 300
    step_width: float = 0.02

    def initialize(self):
        if self.de_model not in DE_PARAMS:
            raise ValueError(f"de_model inconnu : {self.de_model}")
        super().initialize()

    def get_helper_theories(self):
        self._camb_transfers = DETransfers(
            self, "camb.transfers", {"stop_at_error": self.stop_at_error},
            timing=self.timer)
        self._camb_transfers.requires = self._transfer_requires
        return {"camb.transfers": self._camb_transfers}


# Alias : les noms de classe des anciens YAML continuent de fonctionner.
class TachyonDE_Exponential(DarkEnergyCAMB):
    de_model: str = "exponential"


class TachyonDE_InversePower(DarkEnergyCAMB):
    de_model: str = "inverse_power"


class FreeformWz(DarkEnergyCAMB):
    de_model: str = "freeform"


class LCDMReference(DarkEnergyCAMB):
    """Même chemin de code que les autres, sans DE tabulée (w = -1)."""
    de_model: str = "lcdm"
