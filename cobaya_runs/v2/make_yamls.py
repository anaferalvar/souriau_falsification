"""Génère les YAML v2 (après correctif d'injection).

  python cobaya_runs/v2/make_yamls.py --profile full  --packages ~/cobaya_packages
  python cobaya_runs/v2/make_yamls.py --profile local --packages ~/cobaya_packages

profile full  : combinaison canonique de l'amendement-01 (DESI DR2 BAO + Planck 2018
                lite + lowl TT/EE + SNe). À lancer sur serveur (~4-5 s par point).
profile local : DESI DR2 BAO + SNe + prior BBN sur ombh2, sans CMB. Validation rapide
                (~0.05 s par point) : NE remplace PAS le verdict canonique.
"""
import argparse, os, yaml

SNE = {"pantheon": "sn.pantheonplus", "des": "sn.desy5"}   # DES-SN5YR, comme pré-enregistré
MODELS = {
    "lcdm":     ("src.tachyon.cobaya_interface.LCDMReference", {}),
    "expo":     ("src.tachyon.cobaya_interface.TachyonDE_Exponential", {
        "tachyon_lambda":   {"prior": {"min": 0.01, "max": 3.0}, "ref": {"dist": "norm", "loc": 0.5, "scale": 0.1},
                             "proposal": 0.1, "latex": r"\lambda"},
        "tachyon_phi_init": {"prior": {"min": 0.0, "max": 1.5}, "ref": {"dist": "norm", "loc": 0.3, "scale": 0.05},
                             "proposal": 0.05, "latex": r"\phi_i"}}),
    "invpower": ("src.tachyon.cobaya_interface.TachyonDE_InversePower", {
        "tachyon_n":        {"prior": {"min": 0.1, "max": 4.0}, "ref": {"dist": "norm", "loc": 0.5, "scale": 0.1},
                             "proposal": 0.1, "latex": "n"},
        "tachyon_phi_init": {"prior": {"min": 0.1, "max": 2.0}, "ref": {"dist": "norm", "loc": 1.0, "scale": 0.05},
                             "proposal": 0.05, "latex": r"\phi_i"}}),
    "freeform": ("src.tachyon.cobaya_interface.FreeformWz", {
        f"w_bin{i}": {"prior": {"min": -3.0, "max": 1.0}, "ref": {"dist": "norm", "loc": -1.0, "scale": 0.05},
                      "proposal": 0.05, "latex": f"w_{i}"} for i in range(1, 5)}),
}
SEEDS = {"lcdm": 40, "expo": 42, "invpower": 43, "freeform": 44}


def common_params(profile):
    p = {"H0": {"prior": {"min": 60, "max": 75}, "ref": {"dist": "norm", "loc": 68, "scale": 0.5},
                "proposal": 0.3, "latex": "H_0"},
         "omch2": {"prior": {"min": 0.08, "max": 0.16}, "ref": {"dist": "norm", "loc": 0.119, "scale": 0.001},
                   "proposal": 0.001, "latex": r"\Omega_c h^2"},
         "mnu": 0.06,
         "Omega_m": {"derived": "lambda omch2, ombh2, H0: (omch2 + ombh2 + 0.06/93.14) / (H0/100)**2",
                     "latex": r"\Omega_m"},
         "w_today": {"latex": "w_0"}, "phidot2_max": {"latex": r"\dot\phi^2_{max}"}}
    if profile == "local":
        # prior BBN (Schöneberg 2024, utilisé par DESI DR2 « BAO+BBN »)
        p["ombh2"] = {"prior": {"dist": "norm", "loc": 0.02218, "scale": 0.00055},
                      "ref": 0.0222, "proposal": 0.0003, "latex": r"\Omega_b h^2"}
    else:
        p["ombh2"] = {"prior": {"min": 0.018, "max": 0.026}, "ref": {"dist": "norm", "loc": 0.02237, "scale": 0.00015},
                      "proposal": 0.0002, "latex": r"\Omega_b h^2"}
        p["As"] = {"prior": {"min": 5e-10, "max": 1e-8}, "ref": {"dist": "norm", "loc": 2.1e-9, "scale": 2e-11},
                   "proposal": 2e-11, "latex": "A_s"}
        p["ns"] = {"prior": {"min": 0.88, "max": 1.06}, "ref": {"dist": "norm", "loc": 0.9649, "scale": 0.004},
                   "proposal": 0.004, "latex": "n_s"}
        p["tau"] = {"prior": {"min": 0.01, "max": 0.15}, "ref": {"dist": "norm", "loc": 0.055, "scale": 0.006},
                    "proposal": 0.006, "latex": r"\tau"}
    return p


def build(model, sne, profile, packages):
    theory_cls, de_params = MODELS[model]
    lik = {"bao.desi_dr2": None, SNE[sne]: None}
    if profile == "full":
        lik.update({"planck_2018_highl_plik.TTTEEE_lite_native": None,
                    "planck_2018_lowl.TT": None, "planck_2018_lowl.EE": None})
    params = common_params(profile)
    params.update(de_params)
    name = f"{profile}_{model}_{sne}"
    return name, {
        "output": f"cobaya_runs/v2/chains/{name}/{name}",
        "packages_path": packages,
        "theory": {theory_cls: {"stop_at_error": False,
                                "extra_args": {"num_massive_neutrinos": 1, "neutrino_hierarchy": "degenerate"}}},
        "likelihood": lik,
        "params": params,
        "sampler": {"mcmc": {"Rminus1_stop": 0.01, "Rminus1_cl_stop": 0.2, "learn_proposal": True,
                             "max_tries": 10000, "seed": SEEDS[model] + (0 if sne == "pantheon" else 100),
                             "drag": profile == "full", "oversample_power": 0.4}},
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", choices=["full", "local"], required=True)
    ap.add_argument("--packages", required=True, help="chemin ABSOLU des cobaya_packages sur la machine qui lance")
    a = ap.parse_args()
    out = os.path.join(os.path.dirname(__file__), a.profile)
    os.makedirs(out, exist_ok=True)
    for model in MODELS:
        for sne in SNE:
            name, info = build(model, sne, a.profile, a.packages)
            with open(os.path.join(out, name + ".yaml"), "w") as f:
                yaml.safe_dump(info, f, sort_keys=False, allow_unicode=True)
            print("écrit", name)
