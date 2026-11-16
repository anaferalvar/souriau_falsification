"""Figure avant/après du correctif : reconstruction libre w(z). Usage : python -m src.analysis.plot_wz local"""
import sys, json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from src.analysis.verdict import load_chain, wmean
from src.tachyon.potentials import ExponentialPotential
from src.tachyon.rolling_tachyon import build_tachyon_background

profile = sys.argv[1] if len(sys.argv) > 1 else "local"
Z = [0, 0.3, 0.6, 1.0, 1.5]
V1 = {"Pantheon+": ([-0.9946, -1.0158, -0.9528, -0.9444], [1.1164, 1.1402, 1.1439, 1.1742]),
      "DES-SN5YR": ([-1.0749, -0.9566, -1.0073, -1.0206], [1.1722, 1.1587, 1.1431, 1.1663])}
r = json.load(open(f"results/v2/{profile}/verdict.json"))
col = {"Pantheon+": "#D85A30", "DES-SN5YR": "#185FA5"}
off = {"Pantheon+": -0.03, "DES-SN5YR": 0.03}
v2 = {k: (list(v["w_mean"].values()), list(v["w_std"].values())) for k, v in r["sigma_phi"].items()}
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
for a, (title, data) in zip(ax, [("v1 — DE jamais injectée : postérieure = prior (±1,15)", V1),
                                 (f"v2 — après correctif ({profile})", v2)]):
    for k, (m, s) in data.items():
        zc = np.array([(Z[i] + Z[i + 1]) / 2 for i in range(4)]) + off[k]
        a.errorbar(zc, m, yerr=s, fmt="o", color=col[k], capsize=4, label=k)
    a.axhline(-1, color="k", ls=":", lw=1)
    a.axhspan(-4, -1, color="red", alpha=0.05)
    for e in Z[1:-1]:
        a.axvline(e, color="gray", lw=0.5, alpha=0.4)
    a.set_title(title, fontsize=10); a.set_xlabel("z"); a.set_ylim(-2.6, 0.3); a.set_xlim(0, 1.5)
sn_key = {"DES-SN5YR": "des", "Pantheon+": "pantheon"}
for sn, s in sn_key.items():
    lam = r["cells"][f"exponentiel × {sn}"]["tachyon_lambda"][0]
    df = load_chain(f"cobaya_runs/v2/chains/{profile}_expo_{s}/{profile}_expo_{s}")
    bg = build_tachyon_background(ExponentialPotential(1, lam), wmean(df, "Omega_m"), 0.0)
    ax[1].plot(bg["z"], bg["w"], color=col[sn], lw=1.2, alpha=0.8, ls="--")
ax[1].plot([], [], "k--", lw=1, label="tachyon exp. (λ moyen)")
ax[0].set_ylabel("w(z)")
for a in ax:
    a.legend(fontsize=8, loc="lower left")
ax[1].text(0.02, 0.05, "zone rouge : w < −1 (fantôme, interdit au tachyon)", fontsize=7, color="gray",
           transform=ax[1].transAxes, ha="left", va="bottom") if False else None
plt.tight_layout()
plt.savefig(f"results/v2/{profile}/figures/wz_avant_apres.png", dpi=150)
print("ok")
