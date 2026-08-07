for level in 15_20 15_300 300_700 700_1850 1800_1850
do
  declare ingest=$(sbatch --parsable ohc_ingest.slurm $level)
  sbatch --dependency afterok:$ingest verify_store.slurm $level
  declare publish=$(sbatch --parsable --dependency afterok:$ingest publish.slurm $level)
  sbatch --dependency afterok:$publish verify_publish.slurm $level
done
