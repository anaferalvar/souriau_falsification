# Audit du 16/11/2026 — correctifs avant dévoilement

Branche : `fix/de-injection`. Ce document liste ce qui était faux dans la v1, comment
c'est corrigé, et comment le vérifier. Les résultats v1 (`results/sigma_phi.json`,
`results/verdict_*.{csv,json}`) sont **invalides** et conservés uniquement pour la traçabilité.

## B1 — L'énergie noire n'était jamais injectée dans CAMB (critique)

- **Symptôme.** Postérieures des paramètres DE = priors : bins w(z) libres à −1,0 ± 1,14
  (prior plat [−3, 1] : −1 ± 1,155) ; run01 λ = 1,62 ± 0,85 et φ_i = 0,76 ± 0,43
  (priors : 1,5 ± 0,86 et 0,75 ± 0,43) ; corrélation λ–χ² = −0,01 ; ΔDIC ≈ 0.
- **Cause.** `cobaya_interface.py` surchargeait `set_camb_params`, méthode qui n'existe pas
  dans le wrapper CAMB de cobaya (vérifié sur 3.6.2 : jamais appelée). Toutes les chaînes
  « tachyon » et « freeform » étaient des chaînes ΛCDM avec des paramètres muets.
  (Même appelée, la table w(a) étendue au-delà de a = 1 aurait été refusée par CAMB.)
- **Preuve.** χ²(BAO DESI DR2 + Pantheon+) avec le code v1 : identique (16,008 + 1404,759)
  pour λ = 0,3, 1 et 3, alors que λ = 3 donne w₀ = −0,53.
- **Correctif.** Injection dans le helper `CambTransfers` (là où cobaya calcule réellement le
  fond), qui déclare supporter les paramètres DE. Test : `python -m tests.test_injection` —
  écart de χ² de 400 à 880 entre points DE différents ; bins tous à −1 ⇒ χ² ΛCDM exact.

## B2 — Le freeform aurait biaisé le test contre le fantôme

`DarkEnergyFluid` ne peut pas franchir w = −1 ; avec `stop_at_error: False`, les points qui
franchissent auraient reçu une vraisemblance nulle **en silence** — exactement la région
que σ_φ doit sonder. Correctif : `DarkEnergyPPF` pour le freeform ; marches adoucies (tanh,
Δz = 0,02) pour éviter les oscillations de la spline cubique de CAMB.

## B3 — σ_φ ne calculait pas la statistique pré-enregistrée

v1 : seuls les bins sous −1 entraient dans v, χ² à (nb de bins violés) d.d.l., p/2.
Amendement §4 : D = 2(ln L_libre − ln L_{w≥−1}) ⇒ en gaussien, projection de la moyenne
sur le cône w ≥ −1 en métrique C⁻¹ ; loi χ̄² de Chernoff. Correctif :
`src/analysis/verdict.py::sigma_phi`, calibration Monte-Carlo exacte sous w = −1.

## B4 — DIC

v1 : 3⟨χ²⟩ − 2χ²_min, sans les poids. Correct : DIC = 2⟨χ²⟩ − χ²_min, moyennes pondérées.

## B5 — A1 : signe de Ωk inversé et seuil inversé

- DESI DR2 + CMB : Ωk = **+**0,0023 ± 0,0011 (Ωk > 0 = univers **ouvert**). La v1 utilisait
  −0,0024 ± 0,0012 (« negatively curved » lu comme Ωk < 0). Le signal penche ~2σ **contre**
  l'univers fermé de Souriau.
- Le test de falsification calculait Ωk + 3σ > 0 (toujours vrai près de 0) et affichait
  « FALSIFIÉ ». Correct : Ωk − 3σ > 0. Verdict A1 : **NON CONCLUANT** (2,1σ vers ouvert).
- Les autres lignes du forest plot sont à re-vérifier contre leurs sources.

## B6 — Écarts de configuration avec le pré-enregistrement

- Les runs « DES » utilisaient `sn.desdovekie` (recalibration Dovekie 2025), pas
  `sn.desy5` (DES-SN5YR pré-enregistré). v2 : `sn.desy5`. Dovekie = test de robustesse à
  déclarer comme tel.
- Référence ΛCDM : passait par la classe tachyon (λ = 0,001). v2 : `LCDMReference`, même
  chemin de code, sans DE tabulée.
- Seuil de convergence pré-enregistré R−1 < 0,01 ; plusieurs runs v1 étaient entre 0,02 et
  0,04, et les références ΛCDM avaient ~1 000 pas. v2 : R−1 < 0,01 partout.

## B7 — Volet neutrino : non valide

`neutrino_meff.yaml` est un fichier RTF ; `mnu_eff` n'est relié à rien (`mnu` fixé à 0,02) :
le paramètre ne pouvait pas affecter la vraisemblance. Une masse effective négative demande
un code de Boltzmann modifié (Elbers et al. 2025). À refaire ou retirer.

## B8 — Le cross-check GP n'a pas été réalisé

L'amendement §6.3 en fait une condition du verdict global ; le notebook 02 le marque
« reporté ». Tant qu'il manque, le verdict global doit le mentionner.

## Performance

Solveur du fond tachyon : Radau → DOP853 (système non raide), ~16× plus rapide, w(z=0)
identique à 10⁻⁵ près. ~0,05 s par point sans CMB, ~4,5 s avec Planck (2 cœurs).

## Reproduire

```bash
pip install cobaya camb getdist
cobaya-install bao.desi_dr2 sn.pantheonplus sn.desy5 \
  planck_2018_highl_plik.TTTEEE_lite_native planck_2018_lowl.TT planck_2018_lowl.EE -p ~/cobaya_packages
python -m tests.test_injection ~/cobaya_packages        # doit afficher OK x3
python cobaya_runs/v2/make_yamls.py --profile full --packages $HOME/cobaya_packages   # chemin absolu
bash cobaya_runs/v2/run_all.sh full 2 4                 # 2 runs en parallèle, 4 chaînes MPI chacun
python -m src.analysis.verdict full                     # -> results/v2/full/verdict.json
```
