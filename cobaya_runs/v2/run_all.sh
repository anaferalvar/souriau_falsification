#!/usr/bin/env bash
# Lance les 8 runs d'un profil. Usage :
#   bash cobaya_runs/v2/run_all.sh full  [N_PARALLEL] [MPI_PER_RUN]
#   bash cobaya_runs/v2/run_all.sh local [N_PARALLEL]
# Sur serveur, recommandé : MPI_PER_RUN=4 (4 chaînes par run, R-1 inter-chaînes).
# Reprise automatique (--resume) si un run est interrompu.
set -u
PROFILE=${1:?full|local}; NPAR=${2:-2}; NMPI=${3:-1}
cd "$(dirname "$0")/../.."
export PYTHONPATH=$PWD OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
mkdir -p cobaya_runs/v2/logs
ls cobaya_runs/v2/$PROFILE/*.yaml | xargs -P "$NPAR" -I{} bash -c '
  n=$(basename {} .yaml)
  if [ "'"$NMPI"'" -gt 1 ]; then
    mpirun -n '"$NMPI"' cobaya-run --resume {} > cobaya_runs/v2/logs/$n.log 2>&1
  else
    cobaya-run --resume {} > cobaya_runs/v2/logs/$n.log 2>&1
  fi
  echo "[$(date +%H:%M)] terminé : $n (code $?)"'
