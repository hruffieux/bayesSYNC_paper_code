#!/bin/bash

## Misspecification study (Section 5.5): 9 settings x 50 replicates, one (setting, replicate) pair per task, settings varying fastest, so that tasks 1-9 form a complete replicate. --cpus-per-task must match n_cpus in run_misspecification.R.
##
## Written for the CSD3 cluster (Cambridge), where the runs of the paper were
## made: adapt the account, partition and modules to your cluster. Submit from
## the root of the repository, for example
##   sbatch slurm/submit_misspecification.sh
## and set BAYESSYNC_OUTPUT to a folder with enough space.

#SBATCH -J misspec
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH -A mrc-bsu-sl2-cpu
#SBATCH -p icelake
#SBATCH --array=1-450
#SBATCH --output=slurm_misspec_%A_%a.out

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load R/4.3.1-icelake

CMD="Rscript simulations/misspecification/run_misspecification.R ${SLURM_ARRAY_TASK_ID}"

echo -e "JobID: ${SLURM_JOB_ID}, array task: ${SLURM_ARRAY_TASK_ID}\n======"
echo "Time: `date`"
echo "Running on node: `hostname`"
echo "Current directory: `pwd`"
echo -e "\nCPUs per task = ${SLURM_CPUS_PER_TASK}"
echo -e "\nExecuting command:\n==================\n$CMD\n"

eval $CMD
