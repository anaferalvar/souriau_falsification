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
