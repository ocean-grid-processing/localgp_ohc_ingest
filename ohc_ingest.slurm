#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --job-name=ohc_ingest
#SBATCH --output=%x.%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=katie.mills@colorado.edu
#SBATCH --mem=32GB
#SBATCH --exclude=blanca-g6-u21

module purge

cd $SLURM_SUBMIT_DIR
./dist/ohc_ingest config.toml \
  --dir_mean /pl/active/giglio/localGP/OP20260127b/OP20260127b_potentialTemperature_-89.5S_89.5N_20.5W_379.5E_2004_2025_${1}/Results/FullField \
  --dir_ensemble /pl/active/giglio/localGP/OP20260127b/OP20260127b_potentialTemperature_-89.5S_89.5N_20.5W_379.5E_2004_2025_${1}/Results/FullFieldLocalCondSim \
  --dir_out /pl/active/giglio/localGP/OP20260127b/OP20260127b_potentialTemperature_-89.5S_89.5N_20.5W_379.5E_2004_2025_${1}/ohc_ingest \
  --years 2004:2025 --months 1:12 --layer ${1} --tag OP20260112b_repro --no-ensemble
