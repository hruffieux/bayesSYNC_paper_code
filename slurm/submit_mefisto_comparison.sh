#!/bin/bash

## Comparison with MEFISTO (Section 5.4): 2 settings x 5 values of N x 26 seeds, one task each. Serial runs. The memory request of the runs of the paper is not recorded here; MEFISTO with subject-specific trajectories needs more as N grows.
##
## Written for the CSD3 cluster (Cambridge), where the runs of the paper were
## made: adapt the account, partition and modules to your cluster. Submit from
## the root of the repository, for example
##   sbatch slurm/submit_mefisto_comparison.sh
## and set BAYESSYNC_OUTPUT to a folder with enough space.

#SBATCH -J mefisto
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH -A mrc-bsu-sl2-cpu
#SBATCH -p icelake
#SBATCH --array=1-260
#SBATCH --output=slurm_mefisto_%A_%a.out

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load R/4.3.1-icelake

CMD="Rscript simulations/mefisto_comparison/run_mefisto_comparison.R ${SLURM_ARRAY_TASK_ID}"

echo -e "JobID: ${SLURM_JOB_ID}, array task: ${SLURM_ARRAY_TASK_ID}\n======"
echo "Time: `date`"
echo "Running on node: `hostname`"
echo "Current directory: `pwd`"
echo -e "\nCPUs per task = ${SLURM_CPUS_PER_TASK}"
echo -e "\nExecuting command:\n==================\n$CMD\n"

eval $CMD
