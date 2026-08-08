#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --job-name=publish_for_ohc
#SBATCH --output=%x.%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=katie.mills@colorado.edu
#SBATCH --mem=64GB
#SBATCH --exclude=blanca-g6-u21

module purge
module load anaconda/2023.09
conda activate ohc_crosscheck

cd $SLURM_SUBMIT_DIR
python scripts/publish.py \
  /scratch/alpine/wimi7695/gcos_validation/OP20260127b_potentialTemperature_-89.5S_89.5N_20.5W_379.5E_2004_2025_${1}/ohc_gcos_validation_plev${1}.zarr \
  --experiment C --product LocalGP --preset wmo --out /scratch/alpine/wimi7695/gcos_validation/OP20260127b_potentialTemperature_-89.5S_89.5N_20.5W_379.5E_2004_2025_${1}
