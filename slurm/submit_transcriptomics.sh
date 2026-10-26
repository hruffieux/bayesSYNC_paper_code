#!/bin/bash

## Transcriptomic analysis (Section 6.2), on the simulated data unless bool_simulated_data is set to FALSE in the script. Add excl_dead after the seed for the sensitivity analysis without the patients who died. --cpus-per-task must match n_cpus in analysis_transcriptomics.R. The memory request of the run of the paper is not recorded here.
##
## Written for the CSD3 cluster (Cambridge), where the runs of the paper were
## made: adapt the account, partition and modules to your cluster. Submit from
## the root of the repository, for example
##   sbatch slurm/submit_transcriptomics.sh
## and set BAYESSYNC_OUTPUT to a folder with enough space.

#SBATCH -J rnaseq
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --mail-type=END,FAIL
#SBATCH -A mrc-bsu-sl2-cpu
#SBATCH -p icelake
#SBATCH --array=1
#SBATCH --output=slurm_rnaseq_%A_%a.out

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load R/4.3.1-icelake

CMD="Rscript application/transcriptomics/analysis_transcriptomics.R 1"

echo -e "JobID: ${SLURM_JOB_ID}, array task: ${SLURM_ARRAY_TASK_ID}\n======"
echo "Time: `date`"
echo "Running on node: `hostname`"
echo "Current directory: `pwd`"
echo -e "\nCPUs per task = ${SLURM_CPUS_PER_TASK}"
echo -e "\nExecuting command:\n==================\n$CMD\n"

eval $CMD
