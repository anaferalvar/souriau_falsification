# Run 01 — Tachyon exponentiel + DESI DR2 + Pantheon+

## Configuration
- Potentiel : exponentiel V(phi) = V0 * exp(-lambda * phi)
- Données   : BAO DESI DR2 + Pantheon+ (sans CMB Planck)
- Machine   : MacBook Air M1, macOS 12
- Date      : 2026-06-15
- Seed      : 42
- Steps     : 2800 acceptés (convergé R-1 = 0.043)

## Résultats
- H0           = 68.20 ± 0.40
- Omega_m      = 0.303 ± 0.008
- lambda       = 1.62 ± 0.85   (peu contraint sans CMB)
- phi_init     = 0.76 ± 0.43   (peu contraint sans CMB)

## Fichiers
- mac_expo_pantheon.1.txt      : chaîne MCMC
- mac_expo_pantheon.covmat     : matrice de covariance
- mac_expo_pantheon.input.yaml : configuration exacte utilisée

## Limites
Sans CMB Planck, lambda et phi_init restent dégénérés.
Run complet avec Planck prévu sur serveur (run05).
