#!/bin/bash
# psc_schaeffer2020 production run: one 80 GB A100, shared QoS, 48 h (no checkpoints:
# this PSC build has no ADIOS2, so the run must finish inside one allocation).
#
#   sbatch --chdir=<run_dir> psc/job_psc_80gb.sh <params.txt>
#
# <run_dir> must already hold the params file. Outputs: pfd.* (fields), prt.* (every
# 20th particle), ps.*.bin + ps_header.txt (binned phase space, read_ps.py).
#SBATCH -A m5032_g
#SBATCH -C gpu&hbm80g
#SBATCH -q shared
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --gpus-per-task=1
#SBATCH -t 48:00:00
#SBATCH -J psc_sch2020
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
set -euo pipefail
source ~/perlmutter_gpu_warpx.profile >/dev/null 2>&1
PARAMS="${1:?usage: sbatch --chdir=<run_dir> $0 <params.txt>}"
EXE=/pscratch/sd/h/hhelal/psc/build-cuda/src/psc_schaeffer2020
echo "binary: $EXE ($(stat -c %y "$EXE"))"
srun stdbuf -oL -eL "$EXE" "$PARAMS"
