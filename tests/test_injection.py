"""Test de non-régression : le chi2 BAO+SNe doit dépendre des paramètres DE.

Usage : python -m tests.test_injection [packages_path]
"""
import sys, time
from cobaya.model import get_model

PKG = sys.argv[1] if len(sys.argv) > 1 else "~/cobaya_packages"
BASE = {"H0": 68.0, "ombh2": 0.0224, "omch2": 0.119, "mnu": 0.06}


def chi2_for(theory, de_params, points):
    info = {"packages_path": PKG,
            "theory": {theory: {"stop_at_error": True}},
            "likelihood": {"bao.desi_dr2": None, "sn.pantheonplus": None},
            "params": dict(BASE, **{p: {"prior": {"min": lo, "max": hi}}
                                    for p, (lo, hi) in de_params.items()},
                           w_today=None, phidot2_max=None)}
    m = get_model(info)
    out = []
    for pt in points:
        t = time.time()
        r = m.logposterior(pt)
        out.append((pt, round(-2 * sum(r.loglikes), 3), r.derived, round(time.time() - t, 2)))
    return out


if __name__ == "__main__":
    S = "src.tachyon.cobaya_interface."
    cases = [
        (S + "TachyonDE_Exponential", {"tachyon_lambda": (0.01, 3), "tachyon_phi_init": (0, 1.5)},
         [{"tachyon_lambda": l, "tachyon_phi_init": 0.0} for l in (0.3, 1.0, 3.0)]),
        (S + "TachyonDE_InversePower", {"tachyon_n": (0.1, 4), "tachyon_phi_init": (0.1, 2)},
         [{"tachyon_n": n, "tachyon_phi_init": 1.0} for n in (0.5, 2.0, 4.0)]),
        (S + "FreeformWz", {f"w_bin{i}": (-3, 1) for i in range(1, 5)},
         [{"w_bin1": a, "w_bin2": b, "w_bin3": -1.0, "w_bin4": -1.0}
          for a, b in ((-1, -1), (-1.3, -1.3), (-0.7, -0.7))]),
    ]
    ok = True
    for theory, pri, pts in cases:
        res = chi2_for(theory, pri, pts)
        print(theory.split(".")[-1])
        for pt, c2, der, dt in res:
            print(f"   {pt}  chi2={c2}  derived={der}  ({dt}s)")
        spread = max(r[1] for r in res) - min(r[1] for r in res)
        ok &= spread > 1.0
        print(f"   écart de chi2 = {spread:.2f}  ->  {'OK' if spread > 1 else 'ÉCHEC : DE non injectée'}")
    sys.exit(0 if ok else 1)
