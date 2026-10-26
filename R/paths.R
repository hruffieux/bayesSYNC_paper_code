# Shared set-up, sourced at the top of every script. Every script is run from
# the root of the repository, for example
#
#   Rscript simulations/large_scale/run_large_scale.R 1
#
# Output goes to BAYESSYNC_OUTPUT if set, and otherwise to output/ in the
# repository. The simulated application data go to BAYESSYNC_DATA if set, and
# otherwise to data/. Both default folders are ignored by git; the fits of the
# large-scale simulation can reach several tens of Gb per replicate, so on a
# cluster BAYESSYNC_OUTPUT should point to a scratch file system.

if (!file.exists(file.path("R", "paths.R"))) {
  stop("Run the scripts from the root of the repository, ",
       "for example: Rscript simulations/large_scale/run_large_scale.R 1")
}

out_dir <- Sys.getenv("BAYESSYNC_OUTPUT", unset = "output")
data_dir <- Sys.getenv("BAYESSYNC_DATA", unset = "data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# Index of the array task. Taken from the first command-line argument if there
# is one, then from SLURM_ARRAY_TASK_ID, and otherwise set to `default`.
get_job_id <- function(default = 1) {
  job_id <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)[1]))
  if (is.na(job_id)) {
    job_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = default))
  }
  job_id
}
