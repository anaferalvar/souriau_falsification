# Pré-enregistrement — Test du tachyon roulant comme énergie noire
Date de gel: 2026-06-12
Brief de référence : Test du "tachyon roulant" comme solution à l'énergie noire (carnet Souriau v0.12)

---

## Question centrale (irrévocable)

Récemment, les télescopes ont suggéré que ñ'énergie noire s'emballe de manière annormale (franchissant une ligne rouge mathématique appelée "limite fantôme"). Notre question: l'univers fait-il ce truc impossible, ou est-ce que nos équations habituelles (CPL) nous induisent en erreur?

## Modèle testé
Pour vérifier ça, on remplace l'énergie noire standard par un "Champ Tachyonique"

Tachyon roulant de Sen (2002):
  L = -V(phi) * sqrt(1 - phi_dot^2)
  w(a) = phi_dot^2 - 1   =>   w in [-1, 0] TOUJOURS (borne mécanique)

Deux potentiels, testés SÉPARÉMENT:
  - Exponentiel       : V(phi) = V0 * exp(-lambda * phi)
  - Inverse-puissance : V(phi) = V0 * phi^(-n)

## Paramètres libres et priors

| Paramètre        | Prior                  | Note                          |
|------------------|------------------------|-------------------------------|
| H0               | [60, 75]               | standard                      |
| Omega_m          | [0.2, 0.4]             | standard                      |
| Omega_b h^2      | BBN / Planck           | standard                      |
| n_s, A_s, tau    | priors Planck standard | standard                      |
| lambda (expo)    | [0, 3]                 | adimensionné                  |
| n (inv-power)    | [0.1, 4]               | adimensionné                  |
| V0               | dérivé (normalisation) | fixé par Omega_DE aujourd'hui |
| phi_initial      | [0, pi/2] (en unités c=1) | condition initiale du champ |

## Jeux de données — combinaison canonique

Run obligatoire (canonique):
  DESI DR2 BAO + Planck PR4 (CMB + lensing) + Pantheon+
  DESI DR2 BAO + Planck PR4 (CMB + lensing) + DES-SN5YR

RÈGLE : Pantheon+ et DES-SN5YR ne sont JAMAIS combinés dans le même run.
Si les verdicts diffèrent entre les deux, c'est un résultat, pas un bug à corriger.

## Comparaison de modèles

Pour chaque run (tachyon exponentiel, tachyon inverse-puissance):
  - Delta_chi2_MAP = chi2_MAP(tachyon) - chi2_MAP(LCDM)
  - Delta_chi2_MAP vs CPL (w0waCDM) sur les mêmes données
  - DIC : Delta_DIC = DIC(tachyon) - DIC(LCDM)
  - Si possible : facteur de Bayes (nested sampling ou approx. Savage-Dickey)

## Règle de décision (FIGÉE - ne pas modifier après lecture des résultats)

### Critère 1 — Le tachyon reproduit-il les données aussi bien que CPL?
  |Delta_DIC(tachyon vs CPL)| < 5  => "ajustement comparable"

### Critère 2 — Le tachyon a-t-il besoin de w < -1?
  Par construction mécanique, NON (w in [-1,0] toujours).
  On regarde plutôt : phi_dot^2 est-il poussé vers 1 (limite, w->0,
  champ "mort"/poussière) par les données, au point de mal ajuster
  les BAO/SNe à bas z qui préfèrent w proche de -1?

### VERDICT FINAL

- SOUTENUE: Critère 1 vrai (ajustement comparable à CPL ou ΛCDM)
  ET le modèle n'est pas poussé vers w->0 dans la zone z<1
  (c-à-d phi_dot^2 reste << 1 sur tout le range observé).
  => Le franchissement CPL est un artefact ; le tachyon est une
  candidate DE viable. Passage à P1 (neutrinos, déjà greffé en §3bis)
  et préparation Euclid DR1.

- DISFAVORISÉE: Delta_chi2_MAP(tachyon) - Delta_chi2_MAP(CPL) > 10
  (le tachyon ajuste nettement moins bien que CPL)
  OU phi_dot^2 -> 1 est requis par le best-fit sur z<1
  (le champ doit "mourir" pour ajuster, contredisant le comportement
  DE observé à bas z).
  => Cette lecture cosmologique de la strate P^2<0 échoue.
  Le carnet (cœur mathématique) n'est PAS affecté.

- NON CONCLUANT: |Delta_DIC(tachyon vs LCDM)| < 5 ET
  |Delta_DIC(tachyon vs CPL)| < 5, sans signal clair sur phi_dot^2.
  => Attente Euclid DR1 (oct 2026) / DESI DR3.

## Diagnostics supplémentaires obligatoires (pas pour le verdict, pour le rapport)

  - Comparaison Pantheon+ vs DES-SN5YR : verdict identique ? si non,
    le signaler explicitement (cf. brief §3, robustesse).
  - Comparaison potentiel exponentiel vs inverse-puissance : le
    verdict dépend-il de la forme de V(phi) ?
  - Valeur du pivot z_p et w(z_p) pour chaque run, à comparer aux
    valeurs CPL publiées (DESI DR2 II, Table V).

## Hors scope de ce run (rappel brief §5)

  - Axe micro (charge U(1), spin SO(2,1)) : non traité ici.
  - Volet neutrino Sigma_m_nu^eff : traité séparément en
    notebooks/neutrino_meff/, mêmes chaînes, coût marginal.

## Reproductibilité

  - cobaya version : (à renseigner après installation finale)
  - camb version   : (à renseigner après installation finale)
  - seed MCMC      : fixé par run, noté dans chaque .yaml
  - commit de gel de ce document : (à renseigner après git commit)