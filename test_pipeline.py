"""
Test minimal du pipeline cobaya + tachyon.
Lance 10 steps MCMC avec BAO DESI DR2 seulement (sans Planck ni SNe)
pour vérifier que cobaya_interface.py fonctionne avant le vrai run.
"""
import sys, os
sys.path.insert(0, os.path.expanduser("~/Desktop/souriau-falsification"))

from cobaya.run import run

info = {
    "output": "/tmp/test_tachyon_mini",
    "packages_path": os.path.expanduser("~/cobaya_packages"),
    "force": True,

    "sampler": {
        "mcmc": {
            "burn_in": 0,
            "max_tries": 100,
            "Rminus1_stop": 100,   # pas de critère de convergence → juste 10 steps
            "max_samples": 10,
            "seed": 42,
        }
    },

    "theory": {
        "camb": {
            "extra_args": {"num_massive_neutrinos": 1}
        },
        "src.tachyon.cobaya_interface.TachyonDE_Exponential": {
            "stop_at_error": True,
            "Omega_r": 9.15e-5,
            "z_init": 20.0,   # réduit pour le test
            "n_points": 200,
        }
    },

    "likelihood": {
        "bao.desi_dr2": None,   # BAO seulement pour le test rapide
    },

    "params": {
        "H0":    {"prior": {"min": 65, "max": 70}, "ref": 67.5, "proposal": 0.3},
        "ombh2": {"prior": {"min": 0.021, "max": 0.024}, "ref": 0.02237, "proposal": 0.0002},
        "omch2": {"prior": {"min": 0.11, "max": 0.13}, "ref": 0.120, "proposal": 0.001},
        "As":    2.1e-9,
        "ns":    0.965,
        "tau":   0.055,
        "mnu":   0.06,
        "tachyon_lambda": {"prior": {"min": 0.1, "max": 2.0}, "ref": 0.7, "proposal": 0.1},
        "tachyon_phi_init": {"prior": {"min": 0.0, "max": 1.0}, "ref": 0.0, "proposal": 0.05},
        "Omega_m": {
            "derived": "lambda omch2, ombh2, H0: (omch2 + ombh2) / (H0/100)**2",
            "latex": "\\Omega_m"
        },
        "w_today":    {"latex": "w_0"},
        "phidot_max": {"latex": "\\max|\\dot\\phi|"},
    }
}

print("Lancement du test pipeline (10 steps, BAO seulement)...")
updated_info, sampler = run(info)
print("\nPipeline OK !")
print("Paramètres du dernier point :")
chain = sampler.products()["sample"]
print(chain[["H0", "tachyon_lambda", "w_today", "phidot_max"]].tail(3))
