"""
verdict.py — statistiques de décision de l'amendement-01, version corrigée.

Corrections par rapport aux notebooks 02/03 d'origine
-----------------------------------------------------
1. sigma_phi. L'ancien code ne gardait que les bins sous -1 (v_i = max(0, -(w_i+1))),
   formait v^T C^-1 v, puis appliquait un chi2 à (nombre de bins violés) degrés de
   liberté et divisait la p-valeur par 2. Ce n'est pas la statistique pré-enregistrée :
   - D doit être 2(ln L_libre - ln L_{w>=-1}), soit, en approximation gaussienne,
     D = min_{w >= -1} (mu - w)^T C^-1 (mu - w)  (projection sur le cône admissible,
     qui tient compte des corrélations entre bins) ;
   - sa loi sous le null (w = -1, point le moins favorable du bord) est une chi2-barre
     de Chernoff, mélange de chi2 dont les poids dépendent de C. On la calibre ici
     par Monte-Carlo exact : w_hat ~ N(-1, C), on résout le même problème projeté.
   - sigma_phi = Phi^-1(1 - p), unilatéral (pas de p/2).
2. DIC. L'ancien code calculait <chi2> + 2(<chi2> - chi2_min) = 3<chi2> - 2 chi2_min,
   sans les poids de la chaîne. Le DIC est D(theta_hat) + 2 p_D avec
   p_D = <chi2> - D(theta_hat) ; avec theta_hat = meilleur point : DIC = 2<chi2> - chi2_min,
   moyennes pondérées par la colonne `weight` de cobaya.
3. Verdict global. L'ancien test « pathologie » vérifiait la chaîne de caractères
   du dict entier ; il est remplacé par un drapeau booléen par cellule.
"""
import json
import os

import numpy as np
from scipy import stats
from scipy.optimize import lsq_linear

BINS = ["w_bin1", "w_bin2", "w_bin3", "w_bin4"]
BIN_PIVOT = "w_bin2"

SEUIL_SIGMA_HIGH = 3.0
SEUIL_SIGMA_LOW = 2.0
SEUIL_DIC_LOW = 2.0
SEUIL_DIC_HIGH = 5.0
SEUIL_PHIDOT2 = 0.95


# ── chaînes ───────────────────────────────────────────────────────────

def load_chain(root, burn_in=0.3):
    """Charge toutes les chaînes cobaya `root.*.txt` (avec poids), retire le burn-in."""
    import glob
    import pandas as pd
    frames = []
    for f in sorted(glob.glob(root + ".*.txt")):
        if not f.split(".")[-2].isdigit():
            continue
        with open(f) as fh:
            header = fh.readline().lstrip("#").split()
        df = pd.read_csv(f, sep=r"\s+", comment="#", names=header)
        frames.append(df.iloc[int(len(df) * burn_in):])
    if not frames:
        raise FileNotFoundError(root)
    return pd.concat(frames, ignore_index=True)


def wmean(df, col):
    return float(np.average(df[col], weights=df["weight"]))


def wcov(df, cols):
    return np.atleast_2d(np.cov(df[cols].values.T, aweights=df["weight"].values))


# ── DIC ───────────────────────────────────────────────────────────────

def dic(df):
    chi2_mean = wmean(df, "chi2")
    chi2_min = float(df["chi2"].min())
    return {"DIC": 2 * chi2_mean - chi2_min, "chi2_mean": chi2_mean,
            "chi2_min": chi2_min, "p_D": chi2_mean - chi2_min}


# ── sigma_phi ─────────────────────────────────────────────────────────

def _proj_dist(w_hat, A):
    """min_{w >= -1} ||A (w_hat - w)||^2  (A = C^-1/2)."""
    if np.all(w_hat >= -1):
        return 0.0
    res = lsq_linear(A, A @ w_hat, bounds=(-1.0, np.inf))
    r = A @ (w_hat - res.x)
    return float(r @ r)


def sigma_phi(mu, C, n_mc=20000, seed=0):
    """Test joint w_i >= -1 (amendement-01 §4), calibration chi2-barre par Monte-Carlo."""
    mu, C = np.asarray(mu, float), np.asarray(C, float)
    A = np.linalg.cholesky(np.linalg.inv(C)).T          # A^T A = C^-1
    D_obs = _proj_dist(mu, A)
    rng = np.random.default_rng(seed)
    draws = rng.multivariate_normal(-np.ones_like(mu), C, size=n_mc)
    D_null = np.array([_proj_dist(w, A) for w in draws])
    # poids de la chi2-barre : fraction de tirages projetés sur k contraintes actives
    p = float((np.sum(D_null >= D_obs) + 1) / (n_mc + 1)) if D_obs > 0 else 1.0
    s = float(max(0.0, stats.norm.isf(p))) if p < 1 else 0.0
    return {"D": D_obs, "p_value": p, "sigma_phi": s,
            "p_D_null_is_zero": float(np.mean(D_null == 0))}


def freeform_summary(df):
    mu = np.array([wmean(df, b) for b in BINS])
    C = wcov(df, BINS)
    out = sigma_phi(mu, C)
    i = BINS.index(BIN_PIVOT)
    w = df[BIN_PIVOT].values
    out.update({
        "w_mean": dict(zip(BINS, mu.round(4).tolist())),
        "w_std": dict(zip(BINS, np.sqrt(np.diag(C)).round(4).tolist())),
        "P_pivot_below_-1": float(np.sum(df["weight"][w < -1]) / np.sum(df["weight"])),
        "corr": np.round(C / np.sqrt(np.outer(np.diag(C), np.diag(C))), 3).tolist(),
    })
    return out


# ── drapeau P ─────────────────────────────────────────────────────────

def flag_P(df):
    best = df.loc[df["chi2"].idxmin()]
    q95 = float(np.quantile(df["phidot2_max"], 0.95))
    return {"phidot2_max_bestfit": float(best["phidot2_max"]),
            "phidot2_max_q95": q95,
            "flag_P": bool(best["phidot2_max"] > SEUIL_PHIDOT2)}


# ── règles §5 et §6 ───────────────────────────────────────────────────

def verdict_cellule(s_phi, d_dic, P):
    if s_phi >= SEUIL_SIGMA_HIGH:
        return "DISFAVORISÉE", "franchissement physique"
    if P:
        return "DISFAVORISÉE", "pathologie interne"
    if s_phi < SEUIL_SIGMA_LOW and d_dic <= SEUIL_DIC_LOW:
        return "SOUTENUE", "pas de franchissement requis, ajustement ≥ ΛCDM"
    if s_phi < SEUIL_SIGMA_LOW and d_dic >= SEUIL_DIC_HIGH:
        return "DISFAVORISÉE", "forme du modèle"
    return "NON CONCLUANT", "zone intermédiaire"


def verdict_global(cells, sphi_by_sne, gp_discordance=None):
    """cells: {(pot, sne): {'verdict','motif','flag_P','dDIC'}}."""
    s = list(sphi_by_sne.values())
    if all(x >= 3 for x in s):
        return "DISFAVORISÉE (franchissement physique)", "σ_φ ≥ 3 sur les deux catalogues"
    if sum(x >= 3 for x in s) == 1:
        return "NON CONCLUANT", "discordance inter-SNe sur σ_φ"
    if gp_discordance:
        return "NON CONCLUANT", "discordance bins / GP"
    if any(c["flag_P"] for c in cells.values()):
        return "DISFAVORISÉE (pathologie interne)", "au moins une cellule porte P"
    if all(x < 2 for x in s) and all(c["dDIC"] <= SEUIL_DIC_LOW for c in cells.values()):
        return "SOUTENUE", "unanimité des 4 cellules"
    if all(x < 2 for x in s) and all(c["dDIC"] >= SEUIL_DIC_HIGH for c in cells.values()):
        return "DISFAVORISÉE (forme du modèle)", "unanimité des 4 cellules"
    return "NON CONCLUANT", "zone intermédiaire ou discordance entre cellules"


# ── point d'entrée ────────────────────────────────────────────────────

def run(profile, chains_dir="cobaya_runs/v2/chains", out_dir="results/v2"):
    root = lambda m, s: os.path.join(chains_dir, f"{profile}_{m}_{s}", f"{profile}_{m}_{s}")
    sne_names = {"pantheon": "Pantheon+", "des": "DES-SN5YR"}
    pots = {"expo": "exponentiel", "invpower": "inverse-puissance"}
    report = {"profile": profile, "sigma_phi": {}, "cells": {}, "gp_crosscheck": "non réalisé"}

    for s, sname in sne_names.items():
        report["sigma_phi"][sname] = freeform_summary(load_chain(root("freeform", s)))
        report.setdefault("lcdm", {})[sname] = dic(load_chain(root("lcdm", s)))

    cells = {}
    for p, pname in pots.items():
        for s, sname in sne_names.items():
            df = load_chain(root(p, s))
            d = dic(df)
            P = flag_P(df)
            dd = d["DIC"] - report["lcdm"][sname]["DIC"]
            sp = report["sigma_phi"][sname]["sigma_phi"]
            v, m = verdict_cellule(sp, dd, P["flag_P"])
            par = "tachyon_lambda" if p == "expo" else "tachyon_n"
            cells[(pname, sname)] = {
                "sigma_phi": sp, "dDIC": dd, "dchi2_min": d["chi2_min"] - report["lcdm"][sname]["chi2_min"],
                **P, "verdict": v, "motif": m,
                par: [wmean(df, par), float(np.sqrt(wcov(df, [par])[0, 0]))],
                "w_today": [wmean(df, "w_today"), float(np.sqrt(wcov(df, ["w_today"])[0, 0]))],
            }
    g, gm = verdict_global(cells, {k: v["sigma_phi"] for k, v in report["sigma_phi"].items()})
    report["cells"] = {f"{a} × {b}": c for (a, b), c in cells.items()}
    report["verdict_global"] = {"verdict": g, "motif": gm}

    os.makedirs(os.path.join(out_dir, profile), exist_ok=True)
    with open(os.path.join(out_dir, profile, "verdict.json"), "w") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    return report


if __name__ == "__main__":
    import sys
    r = run(sys.argv[1] if len(sys.argv) > 1 else "full")
    print(json.dumps(r, indent=2, ensure_ascii=False))
