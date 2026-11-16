# souriau_falsification

Programme de falsification des conséquences cosmologiques de la lecture « trois strates »
du formalisme de Souriau, confrontées à DESI DR2, Planck 2018 et aux supernovas
(Pantheon+, DES-SN5YR). Règles de décision gelées avant dévoilement :
`notebooks/tachyon_DE/00_preregistration.md` et `preregistration_amendment-01.md`.

> **Audit du 16/11/2026.** La v1 n'injectait jamais l'énergie noire dans CAMB : ses
> chaînes étaient du ΛCDM et son verdict n'a pas de valeur. Détail, correctifs et
> reproduction : **[AUDIT.md](AUDIT.md)**.

## État (v2)

| Volet | Données | Statut |
|---|---|---|
| A1 — Courbure Ωk | DESI DR2 + CMB (littérature) | NON CONCLUANT — Ωk = +0,0023 ± 0,0011, ~2σ vers **ouvert** |
| A2/A5 — Tachyon roulant | BAO + SNe + BBN (validation) | SOUTENUE, voir ci-dessous |
| A2/A5 — Tachyon roulant | + Planck 2018 (canonique) | **À lancer sur serveur** |
| Cross-check GP | — | Non réalisé |
| Volet neutrino | — | Non valide (B7 dans AUDIT.md) |

### Validation sans CMB (`results/v2/local/verdict.json`)

| Cellule | σ_φ | ΔDIC (t − ΛCDM) | Δχ²_min | w₀ | Verdict |
|---|---|---|---|---|---|
| exponentiel × Pantheon+ | 1,15 | −2,9 | −4,9 | −0,91 ± 0,05 | SOUTENUE |
| exponentiel × DES-SN5YR | 0,14 | −10,2 | −12,5 | −0,85 ± 0,05 | SOUTENUE |
| inverse-puissance × Pantheon+ | 1,15 | −3,0 | −4,9 | −0,91 ± 0,04 | SOUTENUE |
| inverse-puissance × DES-SN5YR | 0,14 | −10,2 | −12,3 | −0,85 ± 0,04 | SOUTENUE |

Les données ne réclament pas de franchissement de w = −1, et un champ structurellement
non-fantôme ajuste **mieux** que ΛCDM. Ce n'est **pas** le verdict canonique : la
combinaison pré-enregistrée inclut Planck, qui contraint Ω_m et peut déplacer ces chiffres.

## Structure

```
src/tachyon/        modèle (fond tachyon, potentiels, interface cobaya corrigée)
src/analysis/       σ_φ, DIC, drapeau P, règles §5–§6, figure avant/après
tests/              test de non-régression de l'injection
cobaya_runs/v2/     générateur de YAML (profils full / local), script de lancement
results/v2/         sorties v2  (results/*.json|csv à la racine = v1, invalides)
```
