# Code for "A scalable Bayesian functional factor model for high-dimensional longitudinal molecular data"

Salima Jaoua, Daniel Temko and Hélène Ruffieux.

This repository holds the scripts that produce the numerical results of the paper: the four simulation studies of Section 5 and the two analyses of the COVID-19 cohort in Section 6, with the corresponding material of the Supplementary Material. The method is implemented in the R package [bayesSYNC](https://github.com/hruffieux/bayesSYNC), which these scripts call.

The COVID-19 data cannot be made public. The repository instead provides a script that simulates data emulating them, in the format the analysis scripts read, so that every step of the application can be run.

## Installation

The scripts were run with R 4.3.1 on the cluster and R 4.5 locally. Install the release of bayesSYNC used for the paper, v0.1.0, from GitHub,

```r
# install.packages("remotes")
remotes::install_github("hruffieux/bayesSYNC@v0.1.0")
```

and the following packages from CRAN and Bioconductor.

| Used by | Packages |
|---|---|
| all simulations | `Splinets`, `pROC`, `ellipse`, `gplots`, `matrixStats`, `ggplot2` |
| comparison with MEFISTO | `MOFA2` (Bioconductor) and its Python backend `mofapy2`, `dplyr`, `tidyr`, `knitr`, `kableExtra`, `reshape2` |
| applications | `tidyr`, `dplyr`, `matrixStats`, `gplots`, `ggpubr`, `nlme`, `KEGGREST`, `clusterProfiler`, `fgsea`, `org.Hs.eg.db`, `msigdbr`, `MASS` |

`MOFA2` calls `mofapy2` through `basilisk`, which installs it on first use when `use_basilisk = TRUE`. On the cluster we used `use_basilisk = FALSE`, with `mofapy2` installed in the Python environment seen by R. The comparison with MEFISTO was run with `MOFA2` 1.7.3 and `mofapy2` 0.7.2 (Python 3.12). More recent versions may give different results; in our hands, MEFISTO with subject-specific groups did not run on irregular time grids with `mofapy2` 0.7.3.

## Running the scripts

Every script is run from the root of the repository, for example

```
Rscript simulations/large_scale/run_large_scale.R 151
```

Output goes to the folder named by the environment variable `BAYESSYNC_OUTPUT`, and to `output/` if it is not set. The simulated application data go to `BAYESSYNC_DATA`, by default `data/`. Both default folders are ignored by git.

The scripts that fit the model take the index of an array task as their first argument, so that each replicate of a simulation runs as a separate job. The scripts in `slurm/` submit the arrays used for the paper, on the CSD3 cluster of the University of Cambridge; the account, partition and modules need adapting to other clusters. The large-scale simulation and the comparison with MEFISTO took up to 36 hours per job, and each fit of the large-scale simulation writes an output file of several tens of gigabytes, so these two studies are not meant to be run on a laptop. The same holds for the misspecification study (450 jobs of up to 12 hours on 16 cores) and for the fit of the transcriptomic analysis (16 cores, 256 Gb of memory); the header of each of these scripts gives the resources used. The biomarker panel analysis runs on a laptop.

## Contents

```
R/                       functions shared by the scripts
simulations/
  large_scale/           Sections 5.2 and 5.3
  mefisto_comparison/    Section 5.4
  misspecification/      Section 5.5
application/
  simulate_data/         simulated data emulating the COVID-19 cohort
  biomarker_panel/       Section 6.1
  transcriptomics/       Section 6.2
slurm/                   job arrays for the cluster
```

Each script states at its top what it computes, how to run it and what it writes.

### Large-scale simulation (Sections 5.2 and 5.3)

`run_large_scale.R` fits one replicate of one configuration per array task: N = 100 subjects, p = 20,000 variables, Q = 3 factors, with Q_max = 3 or 5 and vanilla or annealed variational inference, over 50 replicates (tasks 1 to 200). Then

| Script | Produces |
|---|---|
| `plot_figure2_single_replicate.R` | panels of Figure 2, from replicate 1 with Q_max = 5 and annealing |
| `plot_figure_S_loadings.R` | Supplementary Figure S1 |
| `gather_large_scale.R`, then `plot_figure3_annealing.R` | panels of Figure 3 |
| `plot_figure4_elbo.R` | panels of Figure 4 |

### Comparison with MEFISTO (Section 5.4)

`run_mefisto_comparison.R` fits bayesSYNC and one of the two MEFISTO settings to one simulated data set per array task, for N between 20 and 40 (tasks 1 to 260). `summarise_mefisto_comparison.R <N> <setting>` then summarises each value of N and each setting, and writes Figure 5(a) (at N = 30) and the trajectories of Supplementary Material S4. `plot_figure5_table1.R` combines these summaries into Figure 5(b), Supplementary Figure S5 (runtimes on the log scale) and Table 1.

The tables use 25 replicates per value of N, except where a MEFISTO run failed or did not complete within 36 hours. The replicates entering each summary are listed in `summarise_mefisto_comparison.R`, with the reason for each exclusion.

### Misspecification (Section 5.5)

`run_misspecification.R` fits one setting of one replicate per array task: nine settings, a well-specified baseline and two severities of each of four departures from the model, over 50 replicates (tasks 1 to 450). The replicates are paired across settings. `gather_misspecification.R` collects the results and `plot_table2_misspecification.R` writes Table 2 and the residual correlogram described in the text.

### COVID-19 application (Section 6)

First simulate the data,

```
Rscript application/simulate_data/simulate_cohort_data.R
```

which takes a few seconds and writes a biomarker panel of 57 markers on 88 patients and transcriptomic data on 22,354 genes, of which 13,177 are expressed, on 164 patients; 83 patients are in both. The data are generated from the model of the paper, with the population-level parameters estimated on the real data (loadings, mean functions, eigenfunctions, class-specific score distributions and residual variances), which are stored in `application/simulate_data/calibration/`, and with the numbers of patients, samples per patient and sampling days of the real design. `make_calibration.R` documents how these files were made from the fits to the real data; it cannot be run without them. The header of `simulate_cohort_data.R` lists what the simulation does not reproduce. Results on the simulated data illustrate the use of the code and are not expected to reproduce the numbers of the paper.

Biomarker panel (Section 6.1):

| Script | Produces |
|---|---|
| `analysis_biomarker_panel.R` | the fit, the panels of Figures 6 and 7, the scores by severity class and the enrichment analyses on KEGG pathways |
| `analysis_biomarker_panel.R 20 excl_dead` | the sensitivity analysis without the patients who died within the window (Supplementary Material S7) |
| `diagnostics_biomarker_panel.R <run_dir>` | numbers of observations, the association of the scores with severity allowing for them, threshold sensitivity (Section 6.1) and the residual correlogram (Supplementary Material S6) |
| `compare_runs.R <run_dir_a> <run_dir_b>` | correspondence between the factors of two fits, for instance with and without the patients who died |
| `mefisto_base_biomarker_panel.R <run_dir>` | the comparison with MEFISTO-base (Supplementary Material S8) |
| `compare_with_truth_simulated.R <run_dir>` | for a fit to simulated data, the comparison of the estimated loadings, eigenfunctions and scores with the values that generated them, and the figure of Supplementary Material S10 |
| `run_simulated_replicates.sh <n_or_seeds> [n_parallel]` | simulates, fits and compares several data sets, each in its own folder under `output/simulated_replicates/`, with the same seed for the data and the fit; `summarise_simulated_replicates.R` gathers the comparisons |

Supplementary Material S10 shows the simulated data set of seed 8, fitted with seed 8:

```
Rscript application/simulate_data/simulate_cohort_data.R 8
Rscript application/biomarker_panel/analysis_biomarker_panel.R 8
Rscript application/biomarker_panel/compare_with_truth_simulated.R output/<run folder ending in _seed_8>
```

The first command replaces the data simulated with the default seed in `data/`.

Transcriptomics (Section 6.2):

| Script | Produces |
|---|---|
| `analysis_transcriptomics.R` | the fit and its diagnostic figures; add `1 excl_dead` for the sensitivity analysis without the patients who died (Supplementary Material S9) |
| `table3_transcriptomic_factors.R <run_dir>` | Table 3 and the corresponding tables of Supplementary Material S9 |
| `temporal_summaries_transcriptomics.R <run_dir>` | the temporal summaries of Section 6.2 and Supplementary Material S9 |
| `plot_figure_S_transcriptomics_temporal.R <run_dir>` | the figure of Supplementary Material S9 on the temporal structure of the interferon and erythroid factors |
| `plot_figure_S_transcriptomics_factors.R <run_dir>` | the figure of Supplementary Material S9 showing the loadings, eigenfunctions and scores of three factors, the counterpart of Figure 6 |
| `severity_regression_transcriptomics.R <run_dir>` | the proportional-odds regression of severity on the scores of all factors (Section 6.2 and Supplementary Material S9) |
| `cross_analysis_severity.R <panel output.RData> <run_dir>` | agreement between the severity axes of the two analyses (Section 6.2) |
| `replot_transcriptomics_from_slim.R <run_dir>` | the diagnostic figures redrawn from a completed fit, without refitting |

To run the application on the real data, set `bool_simulated_data` to `FALSE` in the two analysis scripts and point `COVID_DATA_DIR` (preprocessed biomarker and clinical data) and `COVID_RNASEQ_FILE` (the transcriptomic EList) to the data. The signs of the factors are not determined by the data; the reversals that give the orientation of Figure 6, and the factor labels of Table 3, refer to the factor numbers of the fits to the real data and are applied to those fits only.

The enrichment analyses download the KEGG pathway and compound lists through the KEGG REST API on first use and cache them in the data folder. Some computing clusters block this API, in which case this step has to be run once on a machine with access.

## Licence

GPL (>= 3), as the bayesSYNC package. See `LICENSE`.
