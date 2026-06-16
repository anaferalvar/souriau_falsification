"""
cobaya_interface.py (v3 — finale)
===================================

Approche correcte pour cobaya + CAMB + w(a) tabulé :

La classe TachyonDE hérite de cobaya.theories.camb.CAMBTheory
(wrapper CAMB de cobaya). On surcharge `set_camb_params` pour
injecter DarkEnergyFluid.set_w_a_table AVANT que CAMB calcule.

C'est la méthode documentée pour les modèles DE personnalisés
dans cobaya (cf. DESI Dark Secrets arXiv:2502.08876).

La logique :
  1. cobaya appelle notre calculate(params)
  2. On calcule w(a) depuis rolling_tachyon.py
  3. On appelle super().set_camb_params(...) pour obtenir
     l'objet CAMBparams standard
  4. On remplace DarkEnergy par DarkEnergyFluid.set_w_a_table(a, w)
  5. CAMB tourne avec ce fond modifié → CMB, BAO, SNe calculés
"""

import numpy as np
import camb
from camb import dark_energy as camb_de

# Import du wrapper CAMB de cobaya
from cobaya.theories.camb import CAMB

from src.tachyon.potentials import ExponentialPotential, InversePowerPotential
from src.tachyon.rolling_tachyon import build_tachyon_background


class _TachyonCAMBBase(CAMB):
    """
    Base : hérite du wrapper CAMB de cobaya et surcharge
    set_camb_params pour injecter w(a) tachyon.
    """

    # Paramètres de configuration (surchargeables dans YAML)
    Omega_r: float = 9.15e-5
    z_init:  float = 50.0
    n_points: int  = 500

    # Nom du paramètre libre tachyon (défini dans les sous-classes)
    _free_param_name: str = ""

    def _get_potential(self, free_param_value, phi_init):
        raise NotImplementedError

    def set_camb_params(self, params_values_dict, **kwargs):
        """
        Surcharge de la méthode cobaya qui prépare les paramètres CAMB.
        On :
          1. Appelle la version parente (ΛCDM standard)
          2. Calcule w(a) depuis le tachyon
          3. Remplace DarkEnergy dans l'objet CAMBparams
        """
        # ── 1. Params CAMB standard (Planck, BBN, etc.) ───────────────
        camb_params = super().set_camb_params(params_values_dict, **kwargs)
        if camb_params is None:
            return None

        # ── 2. Récupérer H0 et Omega_m depuis les params ─────────────
        H0    = params_values_dict.get("H0", 67.5)
        omch2 = params_values_dict.get("omch2", 0.120)
        ombh2 = params_values_dict.get("ombh2", 0.0224)
        h = H0 / 100.0
        Omega_m = (omch2 + ombh2) / h**2

        phi_init       = params_values_dict.get("tachyon_phi_init", 0.0)
        free_param_val = params_values_dict.get(self._free_param_name, 1.0)

        # ── 3. Fond tachyon ───────────────────────────────────────────
        try:
            potential = self._get_potential(free_param_val, phi_init)
            bg = build_tachyon_background(
                potential, Omega_m,
                phi_init=phi_init, x_init=0.0,
                Omega_r=self.Omega_r,
                z_init=self.z_init, n_points=self.n_points
            )
        except Exception as e:
            self.log.debug(f"Fond tachyon échoué : {e}")
            return None

        # Vérification physique
        if np.any(bg["w"] < -1.0 - 1e-4) or np.any(np.abs(bg["phi_dot"]) >= 1.0):
            return None

        # ── 4. Injecter w(a) dans CAMB ───────────────────────────────
        a_arr = bg["a"].copy()
        w_arr = bg["w"].copy()

        # Étendre au-delà de a=1 pour éviter extrapolation
        a_ext = np.concatenate([a_arr, [1.05, 1.5, 2.0]])
        w_ext = np.concatenate([w_arr, [w_arr[-1], w_arr[-1], w_arr[-1]]])

        de = camb_de.DarkEnergyFluid()
        de.set_w_a_table(a_ext, w_ext)
        camb_params.DarkEnergy = de

        # ── 5. Stocker dérivés ────────────────────────────────────────
        self._w_today    = float(w_arr[-1])
        self._phidot_max = float(np.max(np.abs(bg["phi_dot"])))

        return camb_params

    def get_can_provide_params(self):
        base = super().get_can_provide_params() if hasattr(super(), 'get_can_provide_params') else []
        return list(base) + ["w_today", "phidot_max"]

    def get_w_today(self):
        return getattr(self, "_w_today", -1.0)

    def get_phidot_max(self):
        return getattr(self, "_phidot_max", 0.0)


# ─────────────────────────────────────────────────────────────────────────────

class TachyonDE_Exponential(_TachyonCAMBBase):
    """V(phi) = V0 * exp(-lambda * phi)"""

    _free_param_name = "tachyon_lambda"

    params = {
        "tachyon_lambda": {
            "prior": {"min": 0.01, "max": 3.0},
            "ref": 0.7,
            "proposal": 0.1,
            "latex": r"\lambda",
        },
        "tachyon_phi_init": {
            "prior": {"min": 0.0, "max": 1.5},
            "ref": 0.0,
            "proposal": 0.05,
            "latex": r"\phi_i",
        }
    }

    def _get_potential(self, val, phi_init):
        return ExponentialPotential(V0=1.0, lam=val)


class TachyonDE_InversePower(_TachyonCAMBBase):
    """V(phi) = V0 * phi^(-n),  phi > 0"""

    _free_param_name = "tachyon_n"

    params = {
        "tachyon_n": {
            "prior": {"min": 0.1, "max": 4.0},
            "ref": 1.0,
            "proposal": 0.1,
            "latex": r"n",
        },
        "tachyon_phi_init": {
            "prior": {"min": 0.1, "max": 2.0},
            "ref": 1.0,
            "proposal": 0.05,
            "latex": r"\phi_i",
        }
    }

    def _get_potential(self, val, phi_init):
        return InversePowerPotential(V0=1.0, n=val)

class FreeformWz(_TachyonCAMBBase):
    """
    Reconstruction libre (non-CPL) de w(z) en 4 bins indépendants.

    Bornes FIGÉES par l'amendement-01 §3 — NE PAS MODIFIER :
        bin 1 : z in [0.0, 0.3]   (contrainte principale : BGS)
        bin 2 : z in [0.3, 0.6]   (contrainte principale : LRG1/LRG2 — bin PIVOT)
        bin 3 : z in [0.6, 1.0]   (contrainte principale : LRG3+ELG1)
        bin 4 : z in [1.0, 1.5]   (contrainte principale : ELG2, QSO)
        w = -1 fixé pour z > 1.5  (énergie noire sous-dominante)

    Chaque w_i a un prior plat indépendant sur [-3, 1] — aucun lissage
    n'est imposé entre bins (le lissage serait lui-même une hypothèse
    pouvant masquer une excursion réelle sous w=-1).

    Cette classe ne représente PAS un modèle physique — c'est un outil
    statistique pour calculer sigma_phi (amendement §4), indépendant
    de tout choix de potentiel V(phi).
    """

    # Bornes figées — amendement-01 §3
    Z_BINS = [0.0, 0.3, 0.6, 1.0, 1.5]
    Z_MAX_EXTRAPOLATION = 5.0

    params = {
        "w_bin1": {
            "prior": {"min": -3.0, "max": 1.0},
            "ref": {"dist": "norm", "loc": -1.0, "scale": 0.1},
            "proposal": 0.08,
            "latex": r"w_1\,[0,0.3]",
        },
        "w_bin2": {
            "prior": {"min": -3.0, "max": 1.0},
            "ref": {"dist": "norm", "loc": -1.0, "scale": 0.1},
            "proposal": 0.08,
            "latex": r"w_2\,[0.3,0.6]",
        },
        "w_bin3": {
            "prior": {"min": -3.0, "max": 1.0},
            "ref": {"dist": "norm", "loc": -1.0, "scale": 0.1},
            "proposal": 0.08,
            "latex": r"w_3\,[0.6,1.0]",
        },
        "w_bin4": {
            "prior": {"min": -3.0, "max": 1.0},
            "ref": {"dist": "norm", "loc": -1.0, "scale": 0.1},
            "proposal": 0.08,
            "latex": r"w_4\,[1.0,1.5]",
        },
    }

    @staticmethod
    def _build_step_wz(w_vals, z_bins, z_max):
        """
        Construit le tableau (a, w) en escalier pour CAMB.

        w_vals : [w1, w2, w3, w4] sur les bins z_bins.
        Au-delà de z_bins[-1], w est fixé à -1 (amendement §3).
        """
        z_edges_full = list(z_bins) + [z_max]
        w_full = list(w_vals) + [-1.0]

        a_arr, w_arr = [], []
        eps = 1e-8

        # Parcours du z le plus grand au plus petit => a croissant
        for i in reversed(range(len(w_full))):
            z_lo, z_hi = z_edges_full[i], z_edges_full[i + 1]
            a_hi = 1.0 / (1.0 + z_lo)
            a_lo = 1.0 / (1.0 + z_hi)
            a_arr += [a_lo, a_hi - eps]
            w_arr += [w_full[i], w_full[i]]

        a_arr = np.array(a_arr)
        w_arr = np.array(w_arr)
        idx = np.argsort(a_arr)
        a_arr, w_arr = a_arr[idx], w_arr[idx]

        # Extrapolation plate au-delà de a=1
        a_arr = np.append(a_arr, [1.05, 2.0])
        w_arr = np.append(w_arr, [-1.0, -1.0])

        return a_arr, w_arr

    def set_camb_params(self, params_values_dict, **kwargs):
        """
        Surcharge : injecte w(a) en escalier au lieu du fond tachyon
        physique. Pas de résolution d'EDO ici — juste un escalier.
        """
        camb_params = CAMB.set_camb_params(self, params_values_dict, **kwargs)
        if camb_params is None:
            return None

        w_vals = [
            params_values_dict["w_bin1"],
            params_values_dict["w_bin2"],
            params_values_dict["w_bin3"],
            params_values_dict["w_bin4"],
        ]

        a_arr, w_arr = self._build_step_wz(
            w_vals, self.Z_BINS, self.Z_MAX_EXTRAPOLATION
        )

        de = camb_de.DarkEnergyFluid()
        de.set_w_a_table(a_arr, w_arr)
        camb_params.DarkEnergy = de

        # Pas de phi_dot/w_today ici — ce n'est pas un modèle physique
        self._w_today = float(w_vals[0])

        return camb_params

    def get_can_provide_params(self):
        base = super().get_can_provide_params()
        return [p for p in base if p != "phidot_max"]  # pas de phi_dot ici