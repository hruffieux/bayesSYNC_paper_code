#!/bin/bash
## ---------------------------------------------------------------------------
## Replicates of the analysis of the simulated biomarker panel
## (Supplementary Material S10): for each seed, simulates a data set with
## simulate_cohort_data.R, fits it with analysis_biomarker_panel.R and compares
## the fit with the generating values with compare_with_truth_simulated.R.
## The same seed is used to simulate the data and to initialise the fit.
## Each replicate has its own data and output folders under
## output/simulated_replicates/, so the default data/ is left untouched.
## A fit takes about ten minutes on a laptop.
##
## Usage, from the root of the repository:
##   bash application/biomarker_panel/run_simulated_replicates.sh <n_replicates> [n_parallel]
## for seeds 1 to n_replicates, or, for chosen seeds, a comma-separated list
##   bash application/biomarker_panel/run_simulated_replicates.sh 0,99,100 [n_parallel]
## then
##   Rscript application/biomarker_panel/summarise_simulated_replicates.R
## ---------------------------------------------------------------------------

n_repl=${1:?number of replicates, or comma-separated seeds}
n_par=${2:-4}
root=output/simulated_replicates
mkdir -p "$root"

run_one() {
  seed=$1
  dir="$root/seed_$seed"
  mkdir -p "$dir/data" "$dir/output"
  # the KEGG lists are downloaded once into data/ and reused, if present
  cp data/kegg_* "$dir/data/" 2>/dev/null
  export BAYESSYNC_DATA="$dir/data" BAYESSYNC_OUTPUT="$dir/output"
  Rscript application/simulate_data/simulate_cohort_data.R "$seed" > "$dir/log.txt" 2>&1 &&
    Rscript application/biomarker_panel/analysis_biomarker_panel.R "$seed" >> "$dir/log.txt" 2>&1
  # the fit is saved before the figures, so the comparison runs even if a
  # figure of the analysis script fails
  run_dir=$(ls -d "$dir"/output/R1_covid_* 2>/dev/null)
  status=1
  if [ -f "$run_dir/output.RData" ]; then
    Rscript application/biomarker_panel/compare_with_truth_simulated.R "$run_dir" >> "$dir/log.txt" 2>&1
    status=$?
  fi
  # the transcriptomic data are not used here
  rm -f "$dir/data/simulated_transcriptomics.rds"
  echo "seed $seed done (exit $status)"
}
export -f run_one
export root

# seeds whose comparison already exists are skipped, so the script can be rerun
# after a failure
case "$n_repl" in
  *,*) seeds=$(echo "$n_repl" | tr ',' ' ') ;;
  *) seeds=$(seq 1 "$n_repl") ;;
esac
for seed in $seeds; do
  ls "$root"/seed_"$seed"/output/*/compare_with_truth.csv > /dev/null 2>&1 || echo "$seed"
done | xargs -P "$n_par" -I{} bash -c 'run_one {}'
