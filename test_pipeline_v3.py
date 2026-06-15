"""
test_pipeline_v3.py
====================
Test du pipeline avec la bonne architecture : TachyonDE hérite de
cobaya.theories.camb.CAMB et surcharge set_camb_params.
"""
import sys, os
sys.path.insert(0, os.path.expanduser("~/Desktop/souriau-falsification"))

from cobaya.run import run

info = {
    "output": "/tmp/test_tachyon_v3",
    "packages_path": os.path.expanduser("~/cobaya_packages"),
    "force": True,

    "sampler": {
        "mcmc": {
            "burn_in": 0,
            "max_tries": 200,
            "Rminus1_stop": 100,
            "max_samples": 10,
            "seed": 42,
        }
    },

    # On utilise notre classe qui HÉRITE de camb — pas de "theory: camb" séparé
    "theory": {
        "src.tachyon.cobaya_interface.TachyonDE_Exponential": {
            "stop_at_error": True,
            "Omega_r": 9.15e-5,
            "z_init": 20.0,
            "n_points": 200,
            # extra_args passés à CAMB (via la classe parente)
            "extra_args": {
                "num_massive_neutrinos": 1,
            }
        }
    },

    "likelihood": {
        "bao.desi_dr2": None,
    },

    "params": {
        "H0":    {"prior": {"min": 65, "max": 70}, "ref": 67.5, "proposal": 0.3,
                  "latex": "H_0"},
        "ombh2": {"prior": {"min": 0.021, "max": 0.024}, "ref": 0.02237,
                  "proposal": 0.0002, "latex": r"\Omega_b h^2"},
        "omch2": {"prior": {"min": 0.11, "max": 0.13}, "ref": 0.120,
                  "proposal": 0.001, "latex": r"\Omega_c h^2"},
        "As":    2.1e-9,
        "ns":    0.965,
        "tau":   0.055,
        "mnu":   0.06,
        "tachyon_lambda": {
            "prior": {"min": 0.1, "max": 2.0}, "ref": 0.7, "proposal": 0.1,
            "latex": r"\lambda"
        },
        "tachyon_phi_init": {
            "prior": {"min": 0.0, "max": 1.0}, "ref": 0.0, "proposal": 0.05,
            "latex": r"\phi_i"
        },
        "Omega_m": {
            "derived": "lambda omch2, ombh2, H0: (omch2 + ombh2) / (H0/100)**2",
            "latex": r"\Omega_m"
        },
    }
}

print("Test pipeline v3 (hérite de CAMB, set_camb_params)...")
updated_info, sampler = run(info)
print("\nOK !")
products = sampler.products()
chain = products["sample"]
cols = [c for c in ["H0", "Omega_m", "tachyon_lambda", "tachyon_phi_init"]
        if c in chain.columns]
print(chain[cols].tail(3).to_string())