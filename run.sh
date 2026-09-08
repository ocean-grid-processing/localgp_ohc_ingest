for layer in 15_20 15_300 300_700 700_1000 700_1850 1800_1850
do
  meandir=/scratch/alpine/wimi7695/ohc_prod/potential_temperature/OP20260507/${layer}/OP20260507_potentialTemperature_-89.5S_89.5N_20.5W_379.5E_2004_2025_${layer}/Results/FullField
  ensembledir=/scratch/alpine/wimi7695/ohc_prod/potential_temperature/OP20260507/${layer}/OP20260507_potentialTemperature_-89.5S_89.5N_20.5W_379.5E_2004_2025_${layer}/Results/FullFieldLocalCondSim
  outputdir=/scratch/alpine/wimi7695/ohc_prod/results_260824
  tag=260824-OP20260507
  # link to the exact ohc_ingest build you deployed (commit or release URL); stamped as code_version.
  codeversion=https://github.com/argovis/ohc_ingest/commit/REPLACE_WITH_DEPLOYED_SHA
  declare ingest=$(sbatch --parsable ohc_ingest.slurm $meandir $ensembledir $outputdir $layer $tag $codeversion)
  sbatch --dependency afterok:$ingest verify_store.slurm $outputdir $meandir $ensembledir $tag $layer
  declare publish=$(sbatch --parsable --dependency afterok:$ingest publish.slurm $outputdir $tag $layer $codeversion)
  sbatch --dependency afterok:$publish verify_publish.slurm $outputdir/OHC_*${layer}*.nc $meandir $ensembledir
done
