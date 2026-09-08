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
  --dir_mean $1 --dir_ensemble $2 --dir_out $3 --layer ${4} --tag $5 \
  --provenance-link https://github.com/argovis/ohc_ingest/releases/tag/$5 --code-version $6


