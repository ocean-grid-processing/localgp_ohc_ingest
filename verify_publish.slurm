#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --job-name=verify_publish
#SBATCH --output=%x.%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=katie.mills@colorado.edu
#SBATCH --mem=32GB
#SBATCH --exclude=blanca-g6-u21

module purge
module load anaconda/2023.09
conda activate ohc_crosscheck

cd $SLURM_SUBMIT_DIR
python scripts/verify_publish.py --ensemble \
  /scratch/alpine/wimi7695/validation/argo_ohc_Global_2004_2024_${1}/OHC_*.nc \
  /scratch/alpine/wimi7695/validation/argo_ohc_Global_2004_2024_${1}/Results/FullField \
  /scratch/alpine/wimi7695/validation/argo_ohc_Global_2004_2024_${1}/Results/FullFieldLocalCondSim
