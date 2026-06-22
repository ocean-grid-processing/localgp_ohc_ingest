#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --job-name=verify_store
#SBATCH --output=%x.%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=katie.mills@colorado.edu
#SBATCH --mem=32GB
#SBATCH --exclude=blanca-g6-u21

module purge
module load anaconda/2023.09
conda activate ohc_crosscheck

cd $SLURM_SUBMIT_DIR
python scripts/publish.py /scratch/alpine/wimi7695/me4oh/postprocessing/outputs/OP20260507/15_20/ohc_localgp_dev_plev15_20.zarr --experiment C --product LocalGP --preset me4oh  --out /scratch/alpine/wimi7695/me4oh/postprocessing/outputs/OP20260507/15_20/submissions/
