#!/bin/bash

## Large-scale simulation (Sections 5.2 and 5.3): 4 configurations x 50 replicates, one fit per task. --cpus-per-task must match n_cpus in run_large_scale.R.
##
## Written for the CSD3 cluster (Cambridge), where the runs of the paper were
## made: adapt the account, partition and modules to your cluster. Submit from
## the root of the repository, for example
##   sbatch slurm/submit_large_scale.sh
## and set BAYESSYNC_OUTPUT to a folder with enough space.

#SBATCH -J large_scale
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=50
#SBATCH --mem=256G
#SBATCH --mail-type=END,FAIL
#SBATCH -A mrc-bsu-sl2-cpu
#SBATCH -p icelake
#SBATCH --array=1-200
#SBATCH --output=slurm_large_scale_%A_%a.out

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load R/4.3.1-icelake

CMD="Rscript simulations/large_scale/run_large_scale.R ${SLURM_ARRAY_TASK_ID}"

echo -e "JobID: ${SLURM_JOB_ID}, array task: ${SLURM_ARRAY_TASK_ID}\n======"
echo "Time: `date`"
echo "Running on node: `hostname`"
echo "Current directory: `pwd`"
echo -e "\nCPUs per task = ${SLURM_CPUS_PER_TASK}"
echo -e "\nExecuting command:\n==================\n$CMD\n"

eval $CMD
