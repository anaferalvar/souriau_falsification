# Pré-enregistrement — Amendement 01 (autoporté)

## Règle de décision révisée et reconstruction libre de $w(z)$

**Fichier :** `notebooks/tachyon_DE/00_preregistration.md` — amendement 01
**Date :** 2026-06-15 · **Statut :** déposé *avant dévoilement*.

À la date de cet amendement, seuls le test pipeline (10 pas, BAO DESI DR2 seul) et les
validations physiques de fond ont été exécutés. **Aucun des quatre runs MCMC complets
n'a été lancé**, et aucune des statistiques de verdict définies ci-dessous ($\sigma_\phi$,
$\Delta\mathrm{DIC}$, drapeau $P$) n'a été calculée sur la combinaison canonique. Figer ou
modifier ces critères ici **ne viole pas la discipline de gel** : ils restent arrêtés *avant*
toute lecture des résultats qu'ils jugent. À committer **avant** le premier `cobaya-run` de
production.

**Portée.** Cet amendement (i) remplace le bloc « Règle de décision » de
`00_preregistration.md` ; (ii) **fixe entièrement** la reconstruction libre de $w(z)$ et la
statistique $\sigma_\phi$ dont dépend désormais le verdict. Le fichier original n'est pas
modifié ; il porte en tête un renvoi vers le présent amendement. Ce document est **autoporté** :
il contient tout ce qui est nécessaire pour conduire l'analyse et rendre le verdict sans
recours à un autre document.

---

## 1. Contexte figé (restitué pour autonomie — inchangé)

**Question centrale.** Le franchissement $w<-1$ préféré par DESI DR2 est-il *physiquement*
requis, ou un artefact de la rigidité de la paramétrisation CPL ($w_0w_a$) ?

**Modèle physique testé.** Tachyon roulant, $\mathcal{L} = -V(\phi)\sqrt{1-\dot\phi^2}$,
d'où $w = \dot\phi^2 - 1 \in [-1, 0)$ — **structurellement non-fantôme**. Deux potentiels :
exponentiel $V_0 e^{-\lambda\phi}$ et inverse-puissance $V_0\,\phi^{-n}$.

**Données (combinaison canonique, obligatoire, telle que configurée dans les YAML).**
BAO DESI DR2 (`bao.desi_dr2`) + CMB **Planck 2018** (`planck_2018_highl_plik.TTTEEE_lite`
+ `planck_2018_lowl.TT/EE`) + SNe. Deux catalogues SNe traités **séparément** : Pantheon+ et
DES-SN5YR. Quatre runs = 2 potentiels $\times$ 2 catalogues SNe.
*Remarque :* la mention « Planck PR4 » d'une note antérieure était erronée pour ce volet ; le
volet utilise Planck 2018 (PR3). Le lensing PR4 ne concerne que le volet neutrino (§7, hors
périmètre).

**Priors et échantillonnage (inchangés).** 8 paramètres cosmologiques standard + 2 paramètres
tachyon ($\lambda$ ou $n$ ; $\phi_{\rm init}$). Univers plat ($\Omega_K=0$) assumé. Convergence
Gelman–Rubin $R-1 < 0{,}01$. Seeds 42–45.

---

## 2. Motif de la révision (deux défauts de la règle initiale)

La règle gelée initiale disfavorisait sur **« $\Delta\chi^2 > 10$ (tachyon *vs* CPL) »** et
jugeait « soutenue » contre $\Lambda$CDM.

**(a) Circularité.** Prendre CPL comme étalon et disfavoriser le tachyon parce qu'il ajuste
moins bien que CPL, c'est *supposer que la préférence de CPL est réelle*, donc présupposer la
réponse à la question testée. Si CPL est l'artefact (reconstruction libre : l'écart retombe à
$1{,}2$–$1{,}9\sigma$ dans la littérature), « le tachyon perd contre CPL » ne disfavorise rien.

**(b) Co-déclenchement.** Références différentes selon le verdict (« soutenue » *vs* $\Lambda$CDM,
« disfavorisée » *vs* CPL). DESI préférant CPL à $\Lambda$CDM à $\sim 2{,}8$–$4{,}2\sigma$
($\Delta\chi^2 \approx 12$–$20$), si le tachyon — incapable de franchir $-1$ — atterrit près de
$\Lambda$CDM, les deux verdicts se déclenchent *ensemble*.

**Correction.** CPL sort de la logique de verdict (il reste calculé et reporté comme contexte,
jamais comme critère). Le pivot de disfavorisation devient une **reconstruction libre, non-CPL,
de $w(z)$** (§3), résumée en une seule grandeur $\sigma_\phi$ (§4). Toutes les branches partagent
$\sigma_\phi$ comme grandeur primaire et $\Lambda$CDM comme unique référence secondaire.

---

## 3. Reconstruction libre de $w(z)$ — spécification figée

**Estimateur primaire : $w(z)$ en bins.** Choisi plutôt qu'un processus gaussien parce qu'il
n'introduit aucun noyau ni longueur de corrélation à régler — ses degrés de liberté sont tous
explicites.

- **Bins.** $K=4$ intervalles en $z$, bornes **figées** $\{0,\ 0{,}3,\ 0{,}6,\ 1{,}0,\ 1{,}5\}$.
  Au-delà de $z=1{,}5$, l'énergie noire est sous-dominante : on **fixe $w=-1$** pour $z>1{,}5$.
  Ce gel du dernier tronçon évite l'artefact des bins de haut-$z$ quasi non contraints (qui
  peuvent fabriquer *ou* masquer une excursion selon le prior). Les bornes couvrent le domaine
  où BAO DESI (BGS $\sim0{,}3$, LRG $\sim0{,}5$–$0{,}7$, ELG $\sim0{,}8$–$1{,}1$, QSO $\sim1{,}5$)
  et SNe contraignent réellement $w$.
- **Priors.** $w_i$ **indépendants**, plats sur $[-3,\ 1]$ (aucun prior de lissage : le lissage
  est lui-même une hypothèse qui peut atténuer une excursion réelle).
- **Paramètres cosmologiques restants et nuisances** : mêmes priors et mêmes likelihoods que
  les runs tachyon (§1).
- **Implémentation** : $w(z)$ en escalier passé à CAMB via `DarkEnergyFluid.set_w_a_table`
  (même chemin que la classe tachyon), notebook `notebooks/tachyon_DE/02_freeform_wz_reconstruction.ipynb`.
- **Une reconstruction par catalogue SNe** : $\sigma_\phi$ est *model-independent*, il ne dépend
  pas du potentiel ; il a donc **2 valeurs** (Pantheon+, DES-SN5YR).

**Cross-check de robustesse (reporté, hors verdict).** Une reconstruction par **processus
gaussien** sur $w(z)$, noyau exponentiel-carré, fonction moyenne $w=-1$, longueur de corrélation
de prior $\ell_z \sim \mathcal{U}[0{,}2,\ 0{,}8]$, amplitude $\sigma_w \sim \mathcal{U}[0,\ 1]$,
toutes **figées ici**. Sert uniquement à vérifier que le verdict ne tient pas au choix
« bins ». Une **discordance** marquée entre bins et GP sur le franchissement du seuil bascule
le verdict global en *non concluant* (§6).

---

## 4. Grandeurs de décision (définitions figées)

Toutes sur la combinaison canonique, **par catalogue SNe**.

**$\sigma_\phi$ — la donnée exige-t-elle un franchissement passé sous $-1$ ?**
Test **joint** (donc sans cueillette de bin) du vecteur reconstruit $(w_1,\dots,w_4)$ contre la
**contrainte d'inégalité** $w_i \ge -1\ \forall i$ :

$$ D \;=\; 2\big(\ln \mathcal{L}_{\text{libre}} - \ln \mathcal{L}_{\,w \ge -1}\big)\ \ge\ 0, $$

où $\mathcal{L}_{\text{libre}}$ est le maximum de vraisemblance des bins libres et
$\mathcal{L}_{\,w\ge-1}$ celui sous la contrainte. $D>0$ signifie que relâcher le plancher $-1$
améliore l'ajustement, c'est-à-dire que la donnée tire sous $-1$.

*Calibration (primaire, peu coûteuse).* Sous l'approximation gaussienne-linéaire, à partir de la
moyenne et de la covariance $C$ de la postérieure des bins (issues d'un seul MCMC de la
reconstruction libre), $D$ suit une loi **$\bar\chi^2$ de Chernoff** (mélange de $\chi^2$ dont
les poids se calculent depuis $C$ et le cône $w_i\ge-1$) — traitement standard d'un test à
frontière d'inégalité (cf. $\Sigma m_\nu\ge0$, $A_{\rm lens}$). $\sigma_\phi$ = l'équivalent
gaussien unilatéral de la $p$-valeur de $D$ sous cette loi.
*Robustesse (reportée).* Recalibrer $D$ par Monte-Carlo sous le null $w_i=-1$ (réalisations
passées dans le pipeline complet) si la postérieure est nettement non gaussienne.

**Secondaire 1-D reporté (intuition, sans multiplicité)** : $P\!\left(w < -1\right)$ dans le bin
le mieux contraint, **désigné a priori** comme $z\in[0{,}3,\ 0{,}6]$ (LRG).

**$\Delta\mathrm{DIC}_{t/\Lambda} = \mathrm{DIC}(\text{tachyon}) - \mathrm{DIC}(\Lambda\mathrm{CDM})$**,
*positif $\Rightarrow$ $\Lambda$CDM favorisé* ; **une valeur par run** (4). Contrôles reportés :
$\Delta\mathrm{AIC}$, $\Delta\mathrm{BIC}$, $\log B$ si faisable.

**$P$ — pathologie interne (par run)** : la postérieure pousse
$\max_{z\in[0,2]}\dot\phi^2(z) > 0{,}95$ (donc $w\to0$, comportement de poussière, le champ cesse
d'agir comme énergie noire), **ou** échec de *shooting* / non-convergence non résolu.

---

## 5. Règle par cellule (potentiel $\times$ SNe)

Branches **ordonnées** ; la première qui s'applique fixe le verdict (exclusivité mutuelle
garantie).

1. **$\sigma_\phi \ge 3$** → **DISFAVORISÉE — *franchissement physique*.**
   La reconstruction libre (hors CPL) impose $w<-1$ de façon robuste ; un modèle structurellement
   $w\ge-1$, dont le tachyon roulant, ne peut le porter. *Verdict fort visé par le projet.*

2. *sinon si* **$P$** → **DISFAVORISÉE — *pathologie interne*.**
   Le champ dégénère (poussière) ou ne converge pas : la lecture champ ne tient pas comme énergie
   noire, indépendamment de la donnée.

3. *sinon si* **$\sigma_\phi < 2$ et $\Delta\mathrm{DIC}_{t/\Lambda} \le +2$** → **SOUTENUE.**
   Aucun franchissement robuste requis **et** un modèle physique $w\ge-1$ ajuste $\Lambda$CDM à
   l'équivalence ou mieux : signature de l'artefact CPL, dont le tachyon est une réalisation viable.

4. *sinon si* **$\sigma_\phi < 2$ et $\Delta\mathrm{DIC}_{t/\Lambda} \ge +5$** →
   **DISFAVORISÉE — *forme du modèle*.**
   La donnée ne réclame pas de phantom, mais le tachyon perd nettement contre $\Lambda$CDM :
   l'échec est la **forme** de ce $w(z)$ *thawing*, **non** la borne $-1$. À ne *pas* lire comme
   franchissement physique (distinct de la branche 1).

5. *sinon* → **NON CONCLUANT** (en attente d'Euclid DR1 / DESI DR3).
   Couvre $2\le\sigma_\phi<3$ (phantom suggéré, non robuste) et $\sigma_\phi<2$ avec
   $2<\Delta\mathrm{DIC}_{t/\Lambda}<5$ (tachyon faiblement disfavorisé).

---

## 6. Verdict global (concordance sur les 4 cellules)

$\sigma_\phi$ n'ayant que 2 valeurs (une par SNe). Évalué dans l'ordre :

1. $\sigma_\phi \ge 3$ pour **les deux** catalogues SNe → **DISFAVORISÉE (franchissement physique)**.
2. *sinon* $\sigma_\phi \ge 3$ sur **un seul** catalogue → **NON CONCLUANT** ; discordance
   inter-SNe **reportée** comme systématique candidate.
3. *sinon* discordance **bins/GP** sur le franchissement → **NON CONCLUANT** (reportée).
4. *sinon* au moins une cellule porte $P$ → **DISFAVORISÉE (pathologie interne)**.
5. *sinon* $\sigma_\phi<2$ pour les deux SNe **et** $\Delta\mathrm{DIC}_{t/\Lambda}\le+2$ dans les
   **4** cellules → **SOUTENUE**.
6. *sinon* $\sigma_\phi<2$ pour les deux SNe **et** $\Delta\mathrm{DIC}_{t/\Lambda}\ge+5$ dans les
   **4** cellules → **DISFAVORISÉE (forme du modèle)**.
7. *sinon* → **NON CONCLUANT** (toute zone intermédiaire ou discordance inter-potentiel ; reportée).

**Asymétrie assumée.** « Soutenue » et « disfavorisée (forme) » exigent l'**unanimité** des 4
cellules ; toute discordance tombe en *non concluant*. C'est délibéré : on ne déclare la viabilité
(ou un échec de forme net) qu'au robuste à travers les nuisances (potentiel, catalogue SNe). Le
prix est qu'il est facile d'atterrir en *non concluant* — c'est le comportement voulu d'un test
qui refuse de sur-conclure, non un défaut.

---

## 7. Volet neutrino ($\Sigma m_\nu^{\rm eff}$, $A_{\rm lens}$) — hors périmètre

Le volet §3 bis du brief (masse de neutrino effective négative ; diagnostic d'origine par
$A_{\rm lens}$) tourne sur les mêmes chaînes mais **n'est pas régi par cet amendement**. Sa règle
de verdict propre (signal absorbé par $A_{\rm lens}$ → systématique ; reporté sur $w(z)$ → énergie
noire ; résiste aux deux → signal neutrino) est à **auditer séparément** pour les mêmes deux
défauts (référence circulaire, branches non disjointes) avant son propre dévoilement. Il utilise
le lensing CMB **Planck PR4**, distinct du CMB de ce volet (PR3).

---

## 8. Ce qui ne change pas

Question centrale ; lagrangien et les deux potentiels ; combinaison canonique multi-sondes
obligatoire (PR3) et les deux catalogues SNe séparés ; priors cosmologiques et tachyon ;
platitude assumée ; convergence $R-1<0{,}01$ ; seeds 42–45 ; interdiction de cherry-picking de
paramétrisation. CPL **demeure exécuté et reporté**, hors logique de verdict.

---

## 9. Livrables attendus (mise à jour des « prochaines étapes »)

- `02_freeform_wz_reconstruction.ipynb` — reconstruction en bins (primaire) + GP (robustesse),
  calcul de $D$, calibration $\bar\chi^2$, sortie $\sigma_\phi$ par catalogue SNe. **À figer et
  committer avant les runs de verdict.**
- `03_verdict.ipynb` — applique les §5–§6 sur les 4 chaînes + les 2 reconstructions libres ;
  reporte $\sigma_\phi$, $\Delta\mathrm{DIC}_{t/\Lambda}$ (et AIC/BIC/$\log B$), $P$, le contexte
  CPL, et le verdict par cellule puis global.
- Correctif documentaire indépendant : la justification de platitude de `rolling_tachyon.py`
  (« en testant le coasting on a assumé plat ») est un non-sequitur ; remplacer par « platitude
  assumée, $\Omega_K=0$ ; conforté a posteriori par A1 ($\Omega_K=+0{,}0023\pm0{,}0011$) ».
  Sans effet sur la règle.

---

*À committer tel quel. Ligne de commit suggérée :*
`prereg(amend-01): CPL hors verdict (circularite + co-declenchement) ; pivot = reconstruction libre w(z) en bins ; sigma_phi via chi-bar-2 ; branches ordonnees disjointes ; CMB=PR3 ; volet neutrino hors perimetre [pre-unblinding 2026-06-15]`
