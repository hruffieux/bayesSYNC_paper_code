## ---------------------------------------------------------------------------
## Transcriptomic analysis of the COVID-19 cohort, Section 6.2 of the paper
## (Table 3; Supplementary Material S9 is drawn from the output of this script
## by the other scripts of this folder).
##
## Whole-blood bulk RNA-seq (PAXgene) data of the same cohort as the biomarker
## panel. The genes detected in at least 20% of the analysed samples are
## retained (13,177 on the real data), and the model is fitted with Q_max = 10,
## L_max = 5 and K = 4 interior knots.
##
## The real data cannot be made public. By default the script runs on the
## simulated data written by application/simulate_data/simulate_cohort_data.R,
## which have the format of the real EList below; set bool_simulated_data to
## FALSE, COVID_RNASEQ_FILE to the real EList and COVID_DATA_DIR to the folder
## of the preprocessed clinical data to run it on the real data.
##
## Meant for a cluster (slurm/submit_transcriptomics.sh): the fit was given 16
## cores, 256 Gb of memory and up to 36 hours. The other scripts of this folder
## only read its output_slim.RData.
##
## Usage, from the root of the repository:
##   Rscript application/transcriptomics/analysis_transcriptomics.R [seed] [excl_dead]
##
## The analysis of the paper is the default (N = 164 on the real data, seed 1),
## with the patients who died within the window retained. Giving excl_dead as
## second argument excludes them, which is the sensitivity analysis of
## Supplementary Material S9 (N = 149 on the real data).
##
## Writes, in its own folder under out_dir: output_slim.RData (read by the
## other scripts of this folder), output.RData (the full fit, large), out.txt
## and the diagnostic figures of section 12.
##
## Analysis set, mirroring the metabolomics application:
##   - clinical severity classes B to E (asymptomatic class A and healthy
##     controls excluded from the model);
##   - first seven weeks after symptom onset (days_thres = 50, i.e. days 1-49).
##
## Input: Hanson_PAXRNASeq_FilteredDGEList_processed.rds, an limma EList with
##   $E       matrix of 22,354 genes x 610 samples, already normalised
##            (edgeR calcNormFactors), voom-transformed, and corrected for the
##            extraction batch with empiricalBayesLM;
##   $weights voom precision weights (not used here, see note below);
##   $targets sample-level annotation (610 rows);
##   $genes   gene-level annotation (Ensembl ID, HGNC symbol, chromosome,
##            length, description);
##   $design  design matrix used by voom; $span, loess span used by voom.
## ---------------------------------------------------------------------------

rm(list = ls())

source("R/paths.R")

bool_dry_run <- F

args <- commandArgs(trailingOnly = TRUE)
seed <- if (length(args) >= 1) as.integer(args[1]) else 1
set.seed(seed)

# TRUE: drop the patients who died within the window (sensitivity analysis);
# FALSE: keep them (analysis reported in the paper)
excl_dead <- length(args) >= 2 && args[2] == "excl_dead"

bool_simulated_data <- TRUE

require(tidyr)
require(dplyr)
require(matrixStats)
require(bayesSYNC)

source("R/fun_application.R")
source("R/plot_utils.R")


## ---------------------------------------------------------------------------
## 1. Load the RNA-seq data
## ---------------------------------------------------------------------------

if (bool_simulated_data) {
  rna_file <- file.path(data_dir, "simulated_transcriptomics.rds")
} else {
  rna_file <- Sys.getenv("COVID_RNASEQ_FILE") # Hanson_PAXRNASeq_FilteredDGEList_processed.rds
}

elist <- readRDS(rna_file)

stopifnot(all(c("E", "targets", "genes") %in% names(elist)))

E_all <- elist$E            # genes x samples, voom-transformed and batch-corrected
W_all <- elist$weights      # voom precision weights, same dimension as E (see below)
df_targets <- elist$targets # one row per sample
df_genes <- elist$genes     # one row per gene

# the column names of E are the sample identifiers, in the order of $targets
stopifnot(identical(colnames(E_all), as.character(df_targets$sampleID)))
stopifnot(identical(rownames(E_all), as.character(df_genes$GeneID)))
stopifnot(identical(dim(W_all), dim(E_all)))
dimnames(W_all) <- dimnames(E_all)

# NOTE on the voom weights. voom estimates the mean-variance trend of the
# log-cpm values and assigns each observation a precision weight inversely
# proportional to the variance predicted from its expected log-count and library
# size, so that observations based on few reads receive a lower weight.
# bayesSYNC has one noise variance per variable, which captures heteroscedasticity
# *between* genes but not the variation in precision *within* a gene across
# samples. Ignoring the weights therefore treats all observations of a given gene
# as equally precise. This is an approximation, mitigated by the fact that lowly
# expressed genes and low-depth samples were already filtered out upstream. The
# diagnostic in Section 4 quantifies how much the weights vary, and
# max_weight_cv below optionally removes the genes for which they vary most.



## ---------------------------------------------------------------------------
## 2. Sample-level annotation, renamed to the conventions used by
##    fun_application.R and plot_utils.R
## ---------------------------------------------------------------------------

df_ann <- data.frame(
  "sample_id_v2" = as.character(df_targets$sampleID),
  "subject_id" = as.character(df_targets$recordID),
  "days_from_sx_or_swab_imputed" = as.numeric(df_targets$days_since),
  "severity" = as.character(df_targets$severity_group),
  "gender" = as.character(df_targets$gender),
  "age" = as.numeric(df_targets$age),
  "bmi" = as.numeric(df_targets$bmi),
  "hosp_group" = as.character(df_targets$hosp_group),
  "resp" = as.character(df_targets$resp),
  "seqbatch" = as.character(df_targets$seqbatch),
  "extractbatch" = as.character(df_targets$extractbatch),
  "hsCRP_mg_l" = as.numeric(df_targets$hsCRP_mg_l),
  stringsAsFactors = FALSE)
rownames(df_ann) <- df_ann$sample_id_v2

# The healthy controls are excluded from the model but retained in df_comb so
# that the reference interquartile band can be drawn in the trajectory plots.
# They come from two sources: 32 samples carry a recordID in a series of their
# own (CV0900 to CV0949) and 28 carry none, with purely numeric sample
# identifiers. None of them has a time since symptom onset, which is what keeps
# them out of the model.
stopifnot(all(is.na(df_ann$days_from_sx_or_swab_imputed[df_ann$severity == "HC"])))

# every sample without a subject identifier is a healthy control
stopifnot(all(df_ann$severity[is.na(df_ann$subject_id)] == "HC"))

# conversely, every patient sample has a subject identifier and a time
stopifnot(all(!is.na(df_ann$subject_id[df_ann$severity != "HC"])))
stopifnot(all(!is.na(df_ann$days_from_sx_or_swab_imputed[df_ann$severity != "HC"])))

# no control shares a subject identifier with a patient, so appending the
# controls to df_comb cannot corrupt the subject-level matching
stopifnot(length(intersect(df_ann$subject_id[df_ann$severity == "HC"],
                           df_ann$subject_id[df_ann$severity != "HC"])) == 0)

# One inconsistency in the annotation: sample CV0244_D28 is assigned to subject
# CV0262. Rather than guessing which of the two identifiers is correct, the
# sample is excluded. Update this vector if the annotation is corrected.
excluded_samples <- c("CV0244_D28")

bool_mismatch <- !is.na(df_ann$subject_id) &
  substr(df_ann$sample_id_v2, 1, nchar(df_ann$subject_id)) != df_ann$subject_id
if (any(bool_mismatch)) {
  cat("Sample identifier(s) not matching the subject identifier: ",
      paste0(df_ann$sample_id_v2[bool_mismatch], " (subject ",
             df_ann$subject_id[bool_mismatch], ")", collapse = ", "), "\n", sep = "")
  if (!all(df_ann$sample_id_v2[bool_mismatch] %in% excluded_samples)) {
    warning("Some inconsistent samples are not listed in excluded_samples.\n")
  }
}

df_ann <- df_ann[!(df_ann$sample_id_v2 %in% excluded_samples), , drop = FALSE]
E_all <- E_all[, df_ann$sample_id_v2, drop = FALSE]
W_all <- W_all[, df_ann$sample_id_v2, drop = FALSE]
cat("Samples excluded on identifier grounds: ", length(excluded_samples), "\n", sep = "")


## ---------------------------------------------------------------------------
## 3. Analysis options
## ---------------------------------------------------------------------------

days_thres <- 50               # first seven weeks after symptom onset (days 1-49)
tnew <- 0:(days_thres-1)       # timepoints to evaluate the prediction at

min_timepoints <- 1            # 1 keeps subjects with a single sample, as in the metabolomics analysis
# 2 restricts to subjects with repeated measurements

selected_severity_groups <- c("B", "C", "D", "E") # symptomatic PCR-positive patients

# Gene filtering. The primary analysis uses the full panel of 22,354 genes:
# gene_filter = "none". Only the near-zero variance safeguard is applied, which
# matters because bool_scale = TRUE rescales every gene to unit variance, so a
# constant or near-constant gene would be amplified into pure noise. Discarding
# non-informative genes is otherwise the job of the spike-and-slab prior, and
# letting it do that job is the point of the analysis.
#
# The filter is independent of bool_dry_run, which only controls maxit and the
# number of cores.

ind_filter <- 1
gene_filter <- c("none",                    # keep everything passing the safeguards (primary analysis)
                 "random",                  # random subset, for local test runs (see below)
                 "expression",              # largest average expression level
                 "variance",                # largest total variance across the analysed samples
                 "within_subject_variance"  # largest average within-subject variance over time
)[ind_filter]

# Cut-off, applied to every filter other than "none": the genes scoring above the
# gene_quantile quantile of the score distribution are kept, so 0.9 keeps the top
# 10%. A larger value means fewer genes. On the full panel of about 22,350 genes,
# 0.9 leaves roughly 2,200 genes, 0.75 about 5,600, 0.5 about 11,200, and 0.995
# about 110, which is a convenient size for a local test run with "random".
# The quantile is taken over the genes that pass the upstream safeguards, so the
# number retained also depends on bool_protein_coding_only and max_weight_cv.
gene_quantile <- 0.5

# Guidance on gene_filter when the full panel is not used:
#  - "random" is the right choice for test runs. A top-variance subset gives an
#    artificially clean picture and is not comparable to the full run, whereas a
#    random subset preserves the realistic mixture of informative genes and noise.
#  - "expression" is the most defensible criterion if a filtered analysis is to
#    be reported, because it does not preselect on variance, which is precisely
#    the quantity the model decomposes. Note that at a high quantile it selects
#    the most abundant transcripts, which are dominated by short non-coding RNAs:
#    fine for a test run, misleading if reported.
#  - "variance" and "within_subject_variance" preselect on the structure being
#    estimated. "within_subject_variance" targets temporal variability, which is
#    what the factor processes capture, and is therefore the most circular of
#    the four. Report either only as an explicitly labelled sensitivity analysis.

bool_protein_coding_only <- FALSE # TRUE keeps only genes with an HGNC symbol
var_thres <- 1e-8                 # drop genes that are constant over the analysed samples

# Detection filter. This is separate from gene_filter above, and it is the one
# that matters. gene_filter ranks genes on a score and keeps a quantile of them;
# the detection filter asks whether a gene is expressed in whole blood at all.
#
# The panel of 22,354 genes was filtered upstream only for very low counts across
# the whole cohort, so it still contains several thousand genes that are not
# expressed in blood: keratins, muscle and cilia genes, 5S rRNA repeats, snRNAs
# and pseudogenes. Their values are close to the detection floor and carry almost
# no signal, but bool_scale rescales every gene to unit variance, which promotes
# that floor noise to the same footing as real expression. Running on the full
# panel produced three factors out of eight loading exclusively on such genes
# (skin and muscle genes, RNU6 and SNORD small RNAs, RNA5S repeats).
#
# The standard remedy in RNA-seq is edgeR::filterByExpr on raw counts. Those are
# not in this file, which holds voom-transformed and batch-corrected values, so
# the criterion here is applied to E directly: a gene is retained if it exceeds
# detect_thres in at least detect_min_prop of the analysed samples. Expressing it
# as a proportion of samples rather than a count keeps a gene expressed in one
# severity class only, which a whole-cohort threshold would discard.
bool_detection_filter <- TRUE

# Threshold on the voom-transformed values. Setting it to NULL takes the mode of
# the low end of the distribution as the detection floor, which adapts to the
# normalisation rather than hard-coding a value; the diagnostic prints what was
# chosen so that it can be reported.
detect_thres <- NULL
detect_min_prop <- 0.20 # roughly the size of the smallest severity class

# Optional safeguard related to the voom weights being ignored: drop the genes
# whose precision weights vary most across samples, since these are the genes
# for which the constant-precision assumption is least tenable. Inf disables it.
max_weight_cv <- Inf

# Display genes by their HGNC symbol rather than by their Ensembl identifier in
# every output. Genes with no symbol keep their Ensembl identifier, and duplicate
# symbols are made unique. The Ensembl identifiers remain available through
# df_genes_sel and gene_id_of_label, which the enrichment analysis uses.
bool_use_symbols <- TRUE


## ---------------------------------------------------------------------------
## 4. Restrict the samples, then select the genes
## ---------------------------------------------------------------------------

df_ann_sub <- df_ann[df_ann$severity %in% selected_severity_groups &
                       !is.na(df_ann$days_from_sx_or_swab_imputed) &
                       df_ann$days_from_sx_or_swab_imputed < days_thres, , drop = FALSE]

cat("Samples retained (severity ", paste0(selected_severity_groups, collapse = ", "),
    ", days < ", days_thres, "): ", nrow(df_ann_sub), "\n", sep = "")
cat("Subjects retained: ", length(unique(df_ann_sub$subject_id)), "\n", sep = "")

E_sub <- E_all[, df_ann_sub$sample_id_v2, drop = FALSE]
W_sub <- W_all[, df_ann_sub$sample_id_v2, drop = FALSE]

# Diagnostic on the voom weights over the analysed samples. The coefficient of
# variation of the weights within a gene measures how far that gene departs from
# the constant-precision assumption made by the model.
weight_cv <- rowSds(W_sub) / rowMeans(W_sub)
cat("\nVoom weights, coefficient of variation within gene (analysed samples):\n")
print(summary(weight_cv))
cat("Genes with CV > 0.5: ", sum(weight_cv > 0.5, na.rm = TRUE), " / ", length(weight_cv),
    "\n\n", sep = "")


## --- Residual batch structure -----------------------------------------------
## E was corrected for the extraction batch with empiricalBayesLM. That
## correction is effective on the samples it was fitted on: on the 114 control
## samples (HC, or classes A and B at 60 days or more) the extraction batch
## explains less gene variance than permuted labels do, which is the expected
## signature of a removed effect.
##
## It transfers only partially to the samples analysed here, classes B to E
## before day 50, which the fitToSamples argument excludes by construction. On
## those samples the extraction batch explains about four times the null level,
## and this survives adjustment for a spline of days_since and for severity, so
## it is not merely biology confounded with batch. Two contributing factors:
## extraction batch 1 is absent from the fitting samples altogether, yet
## accounts for 61 of the 292 analysed samples, so its coefficient could not be
## estimated; and the remaining batches contributed only 10 to 40 fitting
## samples each, which the empirical Bayes shrinkage of empiricalBayesLM pulls
## towards zero, mechanically under-correcting.
##
## The clean fix is upstream: refit empiricalBayesLM on all samples with the
## biological design passed as retainedCovariates. Failing that, Section 12c
## checks whether the residual structure actually reaches the retained factors,
## and bool_cov in Section 8 provides the within-pipeline sensitivity analysis.
##
## The sequencing batch is perfectly nested within the extraction batch on these
## samples, so it carries no information beyond it and is not adjusted for
## separately.
##
## The proportion of variance explained per gene is compared with the value
## obtained after permuting the labels, which gives the null level (k-1)/(n-1).

r2_by_factor <- function(mat, lab) {
  
  lab <- droplevels(as.factor(lab))
  grand <- rowMeans(mat)
  sst <- rowSums((mat - grand)^2)
  
  ssb <- rep(0, nrow(mat))
  for (l in levels(lab)) {
    idx <- lab == l
    m <- rowMeans(mat[, idx, drop = FALSE])
    ssb <- ssb + sum(idx) * (m - grand)^2
  }
  
  ssb / ifelse(sst == 0, NA, sst)
  
}

bool_batch_diagnostic <- TRUE

if (bool_batch_diagnostic) {
  
  set.seed(seed)
  
  cat("--- Proportion of gene variance explained (analysed samples) ---\n")
  for (nm in c("extractbatch", "seqbatch", "severity", "gender")) {
    lab <- df_ann_sub[, nm]
    r2 <- r2_by_factor(E_sub, lab)
    r2_perm <- r2_by_factor(E_sub, sample(lab))
    cat(sprintf("%-14s k = %d  median R2 = %.4f  (permuted %.4f)  genes with R2 > 0.05: %d (permuted %d)\n",
                nm, nlevels(droplevels(as.factor(lab))),
                median(r2, na.rm = TRUE), median(r2_perm, na.rm = TRUE),
                sum(r2 > 0.05, na.rm = TRUE), sum(r2_perm > 0.05, na.rm = TRUE)))
  }
  
  # the two batch variables are nested, hence the redundancy noted above
  cat("\nExtraction batch by sequencing batch:\n")
  print(table(df_ann_sub$extractbatch, df_ann_sub$seqbatch))
  
  # positive control: on the samples EBLM was fitted on, the batch should explain
  # less variance than permuted labels do
  fit_samples <- df_ann$severity == "HC" |
    (df_ann$severity %in% c("A", "B") & !is.na(df_ann$days_from_sx_or_swab_imputed) &
       df_ann$days_from_sx_or_swab_imputed >= 60)
  
  E_fit <- E_all[, df_ann$sample_id_v2[fit_samples], drop = FALSE]
  r2_fit <- r2_by_factor(E_fit, df_ann$extractbatch[fit_samples])
  r2_fit_perm <- r2_by_factor(E_fit, sample(df_ann$extractbatch[fit_samples]))
  cat(sprintf("\nPositive control on the %d EBLM fitting samples: extractbatch median R2 = %.4f (permuted %.4f)\n",
              sum(fit_samples), median(r2_fit, na.rm = TRUE), median(r2_fit_perm, na.rm = TRUE)))
  
  # batches represented in the fit, and in the analysed samples
  cat("Extraction batches among the fitting samples: ",
      paste0(sort(unique(df_ann$extractbatch[fit_samples])), collapse = ", "), "\n", sep = "")
  cat("Extraction batches among the analysed samples: ",
      paste0(sort(unique(df_ann_sub$extractbatch)), collapse = ", "), "\n", sep = "")
  
  # extraction batches follow recruitment order, so they are associated with time
  # since symptom onset: adjusting for batch may remove genuine temporal signal
  cat("\nExtraction batch by time window:\n")
  print(table(df_ann_sub$extractbatch,
              ifelse(df_ann_sub$days_from_sx_or_swab_imputed < 14, "< 14 d", ">= 14 d")))
  cat("\n")
  
}

# Detection floor of the voom-transformed values. The distribution of E over all
# genes and samples is bimodal: a low mode made of undetected genes sitting at
# the floor, and the bulk of expressed genes above it. The floor is taken as the
# position of that low mode, so the threshold follows the normalisation instead
# of being hard-coded.
find_detection_floor <- function(E_sub) {
  
  v <- as.numeric(E_sub)
  v <- v[is.finite(v)]
  
  d <- density(v, n = 1024)
  
  # the low mode is the highest peak below the overall median
  below <- d$x < median(v)
  if (!any(below)) return(quantile(v, 0.1, na.rm = TRUE))
  
  floor_x <- d$x[below][which.max(d$y[below])]
  
  # the cut-off sits just above that mode, at the first local minimum after it
  after <- d$x > floor_x & d$x < median(v)
  if (!any(after)) return(floor_x)
  
  d$x[after][which.min(d$y[after])]
  
}

# Gene selection uses neither the severity class nor any clinical outcome.
select_genes <- function(E_sub, df_ann_sub, df_genes, weight_cv, gene_filter,
                         gene_quantile, bool_protein_coding_only, var_thres,
                         max_weight_cv, bool_detection_filter, detect_thres,
                         detect_min_prop) {
  
  keep <- rowVars(E_sub) > var_thres
  cat("Genes with non-degenerate variance: ", sum(keep), " / ", nrow(E_sub), "\n", sep = "")
  
  if (bool_detection_filter) {
    
    if (is.null(detect_thres)) {
      detect_thres <- find_detection_floor(E_sub)
      cat("Detection floor taken from the low mode of E: ",
          format(detect_thres, digits = 4), "\n", sep = "")
    }
    
    prop_detected <- rowMeans(E_sub > detect_thres)
    detected <- prop_detected >= detect_min_prop
    
    cat("Genes detected above ", format(detect_thres, digits = 4), " in at least ",
        100 * detect_min_prop, "% of the analysed samples: ", sum(detected), " / ",
        nrow(E_sub), " (", round(100 * mean(detected)), "%)\n", sep = "")
    
    # what the discarded genes look like, so that the filter can be checked
    if (any(!detected)) {
      sym_out <- df_genes$hgnc_symbol[match(rownames(E_sub)[!detected], df_genes$GeneID)]
      sym_out <- sym_out[!is.na(sym_out) & sym_out != ""]
      cat("Examples of discarded genes: ",
          paste0(head(sym_out, 12), collapse = ", "), "\n", sep = "")
      cat("Genes discarded with no HGNC symbol: ",
          sum(!detected) - length(sym_out), "\n", sep = "")
    }
    
    keep <- keep & detected
    
  }
  
  if (bool_protein_coding_only) {
    has_symbol <- !is.na(df_genes$hgnc_symbol) & df_genes$hgnc_symbol != ""
    keep <- keep & has_symbol[match(rownames(E_sub), df_genes$GeneID)]
  }
  
  if (is.finite(max_weight_cv)) {
    keep <- keep & !is.na(weight_cv) & weight_cv <= max_weight_cv
  }
  
  cat("Genes passing all safeguards: ", sum(keep), "\n", sep = "")
  E_keep <- E_sub[keep, , drop = FALSE]
  
  if (gene_filter == "none") {
    return(rownames(E_keep))
  }
  
  stopifnot(gene_quantile > 0 & gene_quantile < 1)
  
  if (gene_filter == "random") {
    
    # test runs: a random subset keeps the realistic mixture of informative genes
    # and noise, unlike a subset selected on variance
    n_target <- max(1, round((1 - gene_quantile) * nrow(E_keep)))
    return(sort(sample(rownames(E_keep), n_target)))
    
  } else if (gene_filter == "variance") {
    
    score <- rowVars(E_keep)
    
  } else if (gene_filter == "within_subject_variance") {
    
    # average, over the subjects with repeated measurements, of the variance of
    # each gene across that subject's time points
    subj <- df_ann_sub$subject_id
    subj_repeated <- names(which(table(subj) >= 2))
    stopifnot(length(subj_repeated) > 0)
    
    mat_var <- sapply(subj_repeated, function(ss) rowVars(E_keep[, subj == ss, drop = FALSE]))
    score <- rowMeans(mat_var)
    
  } else if (gene_filter == "expression") {
    
    score <- rowMeans(E_keep)
    
  } else {
    stop("Unknown gene_filter.")
  }
  
  # cut-off on the score itself, so that the value is interpretable and can be
  # reported. Ties at the cut-off are all retained, hence the possible small
  # discrepancy with the nominal fraction.
  cutoff <- quantile(score, probs = gene_quantile, na.rm = TRUE)
  out <- rownames(E_keep)[!is.na(score) & score > cutoff]
  
  cat("Score cut-off at the ", gene_quantile, " quantile of ", gene_filter, ": ",
      format(cutoff, digits = 4), "\n", sep = "")
  
  sort(out)
  
}

set.seed(seed) # for gene_filter = "random"
selected_genes <- select_genes(E_sub, df_ann_sub, df_genes, weight_cv, gene_filter,
                               gene_quantile, bool_protein_coding_only, var_thres,
                               max_weight_cv, bool_detection_filter, detect_thres,
                               detect_min_prop)

cat("Genes selected (", gene_filter,
    ifelse(gene_filter == "none", "", paste0(", quantile ", gene_quantile)),
    "): ", length(selected_genes), "\n", sep = "")

# gene annotation restricted to the selected genes, kept for the interpretation
# of the loadings and for the enrichment analysis
df_genes_sel <- df_genes[match(selected_genes, df_genes$GeneID), , drop = FALSE]
rownames(df_genes_sel) <- df_genes_sel$GeneID

# Display labels. The plotting functions in plot_utils.R use the variable name
# directly as an axis label, so the only way to show symbols instead of Ensembl
# identifiers is to carry the symbols as variable names throughout. Labels are
# sanitised the same way prepare_data_irregular sanitises variable names, so that
# they survive the bool_cov branch unchanged.
make_gene_labels <- function(gene_ids, df_genes_sel) {
  
  sym <- df_genes_sel[gene_ids, "hgnc_symbol"]
  sym[is.na(sym) | sym == ""] <- NA
  
  lab <- ifelse(is.na(sym), gene_ids, sym)
  lab <- gsub("\\+", "", gsub("\\.", "", gsub(" ", "_", gsub("-", "_", lab))))
  lab <- gsub("\\)", "_", gsub("\\(", "_", gsub("\\:", "_", lab)))
  
  make.unique(lab, sep = "_")
  
}

if (bool_use_symbols) {
  
  gene_labels <- make_gene_labels(selected_genes, df_genes_sel)
  
  n_no_symbol <- sum(is.na(df_genes_sel$hgnc_symbol) | df_genes_sel$hgnc_symbol == "")
  cat("Genes displayed by HGNC symbol: ", length(gene_labels) - n_no_symbol, " / ",
      length(gene_labels), " (", n_no_symbol,
      " have no symbol and keep their Ensembl identifier)\n", sep = "")
  
} else {
  gene_labels <- selected_genes
}

# mapping in both directions, used by the enrichment analysis and by the tables
gene_id_of_label <- setNames(selected_genes, gene_labels)
gene_label_of_id <- setNames(gene_labels, selected_genes)
df_genes_sel$label <- gene_labels

all_var <- gene_labels # variable names used throughout, as in covid_study.R


## ---------------------------------------------------------------------------
## 5. Build df_comb, in the wide format expected by fun_application.R
##    (one row per sample, one column per variable)
##    Healthy controls are kept here for the reference band in the plots.
## ---------------------------------------------------------------------------

df_expr <- as.data.frame(t(E_all[selected_genes, , drop = FALSE]))
stopifnot(identical(rownames(df_expr), rownames(df_ann)))
colnames(df_expr) <- gene_labels

df_comb <- cbind(df_ann, df_expr)


## ---------------------------------------------------------------------------
## 6. Clinical information: vital status and covariates
##    df_info comes from the metabolomics preprocessing and carries the hospital
##    outcome. It is needed for (i) the death-related sensitivity analysis and
##    (ii) the death markers in the trajectory plots. The simulated data carry
##    the same information in elist$clinical.
##    Set bool_clinical to FALSE to run without it.
## ---------------------------------------------------------------------------

bool_clinical <- TRUE

if (bool_clinical && bool_simulated_data) {

  df_info <- elist$clinical

  common_subj <- intersect(unique(df_ann_sub$subject_id), unique(df_info$subject_id))
  df_clin_subj <- df_info[!duplicated(df_info$subject_id),
                          c("subject_id", "hospital_outcome",
                            "hospital_outcome_days_from_start_sx_or_first_pos_swab",
                            "last_observation_days_from_start_sx_or_first_pos_swab"),
                          drop = FALSE]
  idx <- match(df_comb$subject_id, df_clin_subj$subject_id)
  df_comb$hospital_outcome <- df_clin_subj$hospital_outcome[idx]
  df_comb$hospital_outcome_days_from_start_sx_or_first_pos_swab <-
    df_clin_subj$hospital_outcome_days_from_start_sx_or_first_pos_swab[idx]
  df_comb$last_observation_days_from_start_sx_or_first_pos_swab <-
    df_clin_subj$last_observation_days_from_start_sx_or_first_pos_swab[idx]
  subj_unknown <- subj_unknown_analysed <- character(0)

} else if (bool_clinical) {

  data_dir_clin <- Sys.getenv("COVID_DATA_DIR")

  version_info <- c("20201123", "20210218", "20210427", "20210427_fu")[4]
  version_ms <- "20210630"
  version_nmr <- "201120"
  version_glyc <- "010721"
  version_lipid <- "20210614"
  version_ques <- "20210505"
  keep_replicates <- FALSE
  
  load(paste0(data_dir_clin, "/",
              "preprocessed_data/metabolites_ms_v", version_ms,
              "_nmr_v", version_nmr,
              "_glyc_v", version_glyc,
              "_lipids_v", version_lipid,
              "_clin_v", version_info,
              "_ques_v", version_ques,
              ifelse(keep_replicates, "_with_replicates", ""),
              ".RData"))
  
  # CHECK: the RNA-seq subject identifiers are the recordID values (CV0002,
  # CV0003, ...). This assumes df_info$subject_id uses the same convention.
  common_subj <- intersect(unique(df_ann_sub$subject_id), unique(df_info$subject_id))
  cat("Subjects matched to the clinical table: ", length(common_subj), " / ",
      length(unique(df_ann_sub$subject_id)), "\n", sep = "")
  if (length(common_subj) == 0) {
    stop("No subject identifier in common with df_info. Check the identifier convention.")
  }
  
  df_clin_subj <- df_info[!duplicated(df_info$subject_id),
                          c("subject_id", "hospital_outcome",
                            "hospital_outcome_days_from_start_sx_or_first_pos_swab",
                            "last_observation_days_from_start_sx_or_first_pos_swab"),
                          drop = FALSE]
  
  idx <- match(df_comb$subject_id, df_clin_subj$subject_id)
  df_comb$hospital_outcome <- df_clin_subj$hospital_outcome[idx]
  df_comb$hospital_outcome_days_from_start_sx_or_first_pos_swab <-
    df_clin_subj$hospital_outcome_days_from_start_sx_or_first_pos_swab[idx]
  df_comb$last_observation_days_from_start_sx_or_first_pos_swab <-
    df_clin_subj$last_observation_days_from_start_sx_or_first_pos_swab[idx]
  
  # Subjects with RNA-seq data but no clinical record have an unknown vital
  # status. They are treated as censored: not known to have died, and followed up
  # to their last observed sample. This is an assumption, not a fact, so the
  # number concerned is reported and these subjects are flagged so that they can
  # be set aside in the death-related analyses.
  subj_unknown <- setdiff(unique(df_comb$subject_id[!is.na(df_comb$subject_id)]),
                          df_clin_subj$subject_id)
  subj_unknown_analysed <- intersect(subj_unknown, unique(df_ann_sub$subject_id))
  
  if (length(subj_unknown_analysed) > 0) {
    cat("Analysed subjects with no clinical record (vital status unknown): ",
        length(subj_unknown_analysed), " / ", length(unique(df_ann_sub$subject_id)),
        "\n", sep = "")
  }
  
} else {
  
  # placeholders so that the plotting functions run: every subject is treated as
  # alive and followed beyond the end of the analysis window
  df_comb$hospital_outcome <- "alive"
  df_comb$hospital_outcome_days_from_start_sx_or_first_pos_swab <- NA
  df_comb$last_observation_days_from_start_sx_or_first_pos_swab <- days_thres + 1
  subj_unknown <- subj_unknown_analysed <- character(0)
  
}

# healthy controls have no time of last observation, and their vital status is
# irrelevant here since they are excluded from the model; the plotting functions
# require a single non-missing value per subject
df_comb$hospital_outcome[df_comb$severity == "HC"] <- "alive"
df_comb$last_observation_days_from_start_sx_or_first_pos_swab[df_comb$severity == "HC"] <- days_thres + 1

# plot_trajectory and plot_factor_trajectory build, for each subject, a status
# and a follow-up time, and require exactly one non-missing value of each:
#   stopifnot(length(status) == 1 & length(time) == 1)
# A missing hospital_outcome makes the status NA, and a missing follow-up time
# leaves the time vector empty, which is what triggers that error. Subject CV0001
# is one such case: recorded as discharged, but with both date fields missing.
# The blocks below fill the gaps so that the failure cannot happen deep inside a
# plotting call, and report what was filled.

# missing outcome: treated as not known to have died
bool_no_outcome <- is.na(df_comb$hospital_outcome)
if (any(bool_no_outcome)) {
  cat("Samples with missing hospital outcome, treated as censored: ",
      sum(bool_no_outcome), " (", length(unique(df_comb$subject_id[bool_no_outcome])),
      " subjects)\n", sep = "")
  df_comb$hospital_outcome[bool_no_outcome] <- "alive"
}

# missing follow-up time for a subject not recorded as dead: censored at the last
# observed sample of that subject, which is known from the RNA-seq data itself
bool_no_time <- !(df_comb$hospital_outcome %in% "dead") &
  is.na(df_comb$last_observation_days_from_start_sx_or_first_pos_swab) &
  !is.na(df_comb$subject_id)

if (any(bool_no_time)) {
  
  last_obs <- tapply(df_comb$days_from_sx_or_swab_imputed[!is.na(df_comb$subject_id)],
                     df_comb$subject_id[!is.na(df_comb$subject_id)],
                     function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
  
  cat("Subjects with no follow-up time, censored at their last observed sample: ",
      length(unique(df_comb$subject_id[bool_no_time])), "\n", sep = "")
  
  df_comb$last_observation_days_from_start_sx_or_first_pos_swab[bool_no_time] <-
    last_obs[df_comb$subject_id[bool_no_time]]
  
}

# a subject recorded as dead but with no date of death cannot be placed in time;
# the plotting functions would fail on it, so this is reported rather than filled
bool_dead_no_date <- df_comb$hospital_outcome %in% "dead" &
  is.na(df_comb$hospital_outcome_days_from_start_sx_or_first_pos_swab)
if (any(bool_dead_no_date)) {
  warning(paste0(length(unique(df_comb$subject_id[bool_dead_no_date])),
                 " subject(s) recorded as dead with no date of death.\n"))
}

# final check, so that any remaining gap is caught here rather than inside a plot
df_check <- df_comb[!is.na(df_comb$subject_id), ]
df_check$time <- ifelse(df_check$hospital_outcome %in% "dead",
                        df_check$hospital_outcome_days_from_start_sx_or_first_pos_swab,
                        df_check$last_observation_days_from_start_sx_or_first_pos_swab)
n_time <- tapply(df_check$time, df_check$subject_id, function(x) length(unique(x[!is.na(x)])))
if (any(n_time != 1)) {
  warning(paste0("Subject(s) without exactly one follow-up time: ",
                 paste0(names(n_time)[n_time != 1], collapse = ", "),
                 ". The trajectory plots will fail on these.\n"))
}
rm(df_check)

## ---------------------------------------------------------------------------
## 7. Severity classes, time window, and the sensitivity analysis excluding
##    subjects who died
## ---------------------------------------------------------------------------

df_comb_no_HC <- choose_severity_groups(df_comb, selected_severity_groups)
df_comb_no_HC <- df_comb_no_HC[df_comb_no_HC$days_from_sx_or_swab_imputed < days_thres, ]

# excl_dead is set from the command line at the top of the script

# "in_window": drop only the subjects who died within the days_thres window, i.e.
#              those whose observed trajectory is truncated by death (default).
# "any":       drop every subject with hospital outcome "dead", including those
#              who died after the end of the window.
dead_definition <- c("in_window", "any")[1]

# Which values of hospital_outcome count as a death. The recorded categories are
# dead, discharged, hospital, other_facility, other_hospital, palliative, and a
# sizeable number of missing values.
#   - "dead" is unambiguous.
#   - "palliative" is not a recorded death but is not a survival either: the
#     subject may have died without the field being updated. Add it here to test
#     the sensitivity of the conclusions to that choice.
#   - a missing outcome is an unknown status, not a survival. Such subjects are
#     counted separately and the death counts are lower bounds.
dead_labels <- c("dead")

if (bool_clinical) {
  
  df_status <- df_clin_subj
  df_status$dead <- df_status$hospital_outcome %in% dead_labels
  df_status$unknown_status <- is.na(df_status$hospital_outcome)
  
  # day of death, counted from symptom onset or first positive swab
  df_status$time_death <- rep(NA_real_, nrow(df_status))
  df_status$time_death[df_status$dead] <-
    df_status$hospital_outcome_days_from_start_sx_or_first_pos_swab[df_status$dead]
  
  # time of the event, used by the plotting functions to cut the fitted curve
  df_status$time_event <- df_status$last_observation_days_from_start_sx_or_first_pos_swab
  df_status$time_event[df_status$dead] <- df_status$time_death[df_status$dead]
  
  # a subject who died after the window has a fully observed trajectory over the
  # window and is not concerned by the informative truncation issue
  df_status$dead_in_window <- df_status$dead & !is.na(df_status$time_death) &
    df_status$time_death < days_thres
  
  if (any(df_status$dead & is.na(df_status$time_death))) {
    warning(paste0(sum(df_status$dead & is.na(df_status$time_death)),
                   " subject(s) recorded as dead but with no date of death; ",
                   "counted as dying outside the window.\n"))
  }
  
  if (excl_dead) {
    
    stopifnot(dead_definition %in% c("in_window", "any"))
    subj_to_drop <- df_status$subject_id[if (dead_definition == "any") df_status$dead else df_status$dead_in_window]
    
    n_before <- length(unique(df_comb_no_HC$subject_id))
    df_comb_no_HC <- df_comb_no_HC[!(df_comb_no_HC$subject_id %in% subj_to_drop), , drop = FALSE]
    n_after <- length(unique(df_comb_no_HC$subject_id))
    
    cat(paste0("Excluding subjects who died (definition: ", dead_definition, "): ",
               n_before - n_after, " subjects removed before data preparation.\n"))
    
  }
  
} else {
  stopifnot(!excl_dead) # the vital status is only available through df_info
  df_status <- NULL
}

# Optional exclusion of named subjects, to check their influence on the fit. In
# the Q = 10, K = 4 fit, CV0047 has scores 13 to 16 standard deviations from
# the other subjects on the leading components of factors 2, 5 and 8, and its
# fitted trajectories reach large values at day 0 although its samples span
# days 21 to 49. The choice is recorded in the name of res_dir. The analysis
# of the paper excludes it; the identifier belongs to the real data.
bool_excl_subjects <- !bool_simulated_data
vec_excl_subjects <- c("CV0047") # only used when bool_excl_subjects is TRUE

if (bool_excl_subjects) {
  stopifnot(length(vec_excl_subjects) > 0)
  missing_subj <- setdiff(vec_excl_subjects, df_comb_no_HC$subject_id)
  if (length(missing_subj) > 0) {
    stop(paste0("Subject(s) to exclude not in the analysis set: ",
                paste(missing_subj, collapse = ", ")))
  }
  df_comb_no_HC <- df_comb_no_HC[!(df_comb_no_HC$subject_id %in% vec_excl_subjects), , drop = FALSE]
  cat(paste0("Excluding named subjects: ", paste(vec_excl_subjects, collapse = ", "), "\n"))
}


## ---------------------------------------------------------------------------
## 8. Optional adjustment for covariates
##    E is already corrected for the extraction batch (empiricalBayesLM), but the
##    correction was fitted on the control samples only and leaves residual batch
##    structure on the samples analysed here, as quantified in Section 4.
##    Setting bool_cov to TRUE regresses the covariates out with a subject-level
##    random intercept, as in the metabolomics analysis.
##
##    Two points to keep in mind before switching this on.
##    - The sequencing batch is perfectly nested within the extraction batch on
##      these samples, so adjusting for both is redundant. Only the extraction
##      batch is used.
##    - The extraction batches follow recruitment order and are therefore
##      associated with time since symptom onset. Adjusting for them can remove
##      genuine temporal signal, which is the object of the analysis. Treat this
##      as a sensitivity analysis rather than as the primary specification, and
##      read it alongside the post-hoc diagnostic in Section 12c.
## ---------------------------------------------------------------------------

bool_cov <- F

adjust_for <- c("extractbatch", "age", "gender") # subset of covariates to regress out

if (bool_cov) {
  
  df_clin <- df_ann[df_ann$sample_id_v2 %in% df_comb_no_HC$sample_id_v2, , drop = FALSE]
  
  covariates <- NULL
  
  if ("gender" %in% adjust_for) {
    df_clin$gender <- as.numeric(factor(df_clin$gender, levels = c("F", "M"))) - 1 # 1 = M, 0 = F
    covariates <- c(covariates, "gender")
  }
  
  if ("age" %in% adjust_for) {
    covariates <- c(covariates, "age")
  }
  
  if ("extractbatch" %in% adjust_for) {
    # categorical: dummy coding, one column dropped to avoid collinearity with
    # the intercept
    mat_batch <- mixOmics::unmap(df_clin$extractbatch)
    colnames(mat_batch) <- paste0("extractbatch_", sort(unique(df_clin$extractbatch)))
    df_clin <- cbind(df_clin, mat_batch[, -1, drop = FALSE])
    covariates <- c(covariates, colnames(mat_batch)[-1])
  }
  
  stopifnot(!is.null(covariates))
  df_clin <- df_clin[, c("sample_id_v2", "subject_id", covariates), drop = FALSE]
  
} else {
  df_clin <- NULL
}


## ---------------------------------------------------------------------------
## 9. Prepare the data for bayesSYNC
## ---------------------------------------------------------------------------

bool_rescale_time <- TRUE
data <- prepare_data_irregular(df_comb_no_HC, all_var, min_nb_timepoints = min_timepoints,
                               days_thres = days_thres, bool_rescale_time = bool_rescale_time,
                               Z = df_clin)

list_time_obs_Y <- enforce_common_grid(data$time_obs, data$Y, min_nb_timepoints = min_timepoints)

time_obs <- list_time_obs_Y$time_obs
Y <- list_time_obs_Y$Y
N <- length(Y)
p <- data$p
N
p

subject_ids <- names(Y)

n_obs <- sapply(time_obs, length)
N_min <- min(n_obs)
N_max <- max(n_obs)
N_min
N_max


## ---------------------------------------------------------------------------
## 10. Description of the analysis set
##     (number of subjects and number of observations per subject)
## ---------------------------------------------------------------------------

# cat("\n=== ANALYSIS SET ===\n")
# cat("Subjects: ", N, "\n", sep = "")
# cat("Genes: ", p, "\n", sep = "")
# cat("Total observations: ", sum(n_obs), "\n", sep = "")
# cat("Observations per subject: min ", N_min, ", median ", median(n_obs),
#     ", max ", N_max, "\n", sep = "")
# cat("\nDistribution of the number of observations per subject:\n")
# print(table(n_obs))
# cat("\nSeverity class:\n")
# print(table(df_comb_no_HC$severity[match(names(Y), df_comb_no_HC$subject_id)]))

if (bool_clinical) {
  
  df_status_Y <- df_status[match(names(Y), df_status$subject_id), , drop = FALSE]
  n_unknown <- sum(names(Y) %in% subj_unknown_analysed)
  
  cat("\nDeath defined as hospital_outcome in: ", paste0(dead_labels, collapse = ", "), "\n", sep = "")
  cat("Deaths (any time): ", sum(df_status_Y$dead, na.rm = TRUE), " / ", N, "\n", sep = "")
  cat("Deaths within the ", days_thres, "-day window: ",
      sum(df_status_Y$dead_in_window, na.rm = TRUE), " / ", N, "\n", sep = "")
  cat("Vital status unknown: ",
      sum(df_status_Y$unknown_status, na.rm = TRUE) + n_unknown, " / ", N,
      " (", sum(df_status_Y$unknown_status, na.rm = TRUE), " with a missing outcome, ",
      n_unknown, " with no clinical record). The two counts above are lower bounds.\n", sep = "")
  
  cat("\nRecorded hospital outcome:\n")
  print(table(df_status_Y$hospital_outcome, useNA = "ifany"))
  
}
cat("====================\n\n")


## ---------------------------------------------------------------------------
## 11. Model settings and inference
## ---------------------------------------------------------------------------

# Number of FPCA components. The per-factor cumulated PVE reaches 99.9 to 100 per
# cent at l = 4 for all five factors, so L = 4 loses essentially nothing and cuts
# the runtime, the dominant cost being O(N * Q * L * K * p) per iteration. Unlike
# the choice of K below, this is a modelling decision and should be stated if it
# differs from the simulations.
L <- 5
Q <- 10

# Number of spline basis functions. Leaving K at NULL applies the package rule
#   K <- max(round(min(median(sapply(time_obs, length) / 4), 40)), 10),
# which is the rule of Ruppert (2002) used in the paper. That rule is calibrated
# on the number of observations along a single curve; here the median number of
# observations per subject is 2, so median(2/4) = 0.5 rounds to 0 and the rule
# bottoms out at its floor, giving K = 10, that is K' = 8 penalised basis
# functions. Combined with the deliberately diffuse half-Cauchy(A = 1e5) prior on
# the smoothing variance, nothing then prevents the eigenfunctions from
# oscillating several times over the 50-day window, which is what was observed.
# Setting K explicitly is the direct remedy, and it also cuts the runtime.
K <- 4

bool_anneal <- T
if (bool_anneal) {
  anneal <- c(1, 1.9, 100)
} else {
  anneal <- NULL
}

# Sparsity of the loadings. Two knobs, of very different standing.
#
# bool_var_spec_prob gives each gene-factor pair its own inclusion probability
# omega_jq instead of one omega_q shared by all genes of a factor. In the shared
# case the update is
#     c_1 = c_0 + sum_j gamma_jq,   d_1 = d_0 + p - sum_j gamma_jq,
# so with about 1800 genes selected out of 4471 the inclusion log-odds threshold
# sits near log(7152) - log(1790) = 1.4. In the variable-specific case each
# omega_jq is updated by a single gamma_jq, giving c_1 close to 1 and d_1 close
# to d_0 = p, hence a threshold near log(p) = 8.4. Switching it on therefore
# tightens the selection by roughly 7.6 log-odds units, which is about what
# d_0 = p^2 would do, but through a documented option with a clear reading: no
# pooling towards a common sparsity level.
bool_var_spec_prob <- FALSE

# d_0_power raises the second parameter of the Beta prior to d_0 = p^d_0_power.
# The package default is 1, that is d_0 = p, so the prior mean inclusion
# probability is 1 / (p + 1), about one active gene per factor a priori: the
# multiplicity adjustment is already there. Raising it should be reported as a
# sensitivity analysis, not used as a tuning knob, since d_0 = p^2 states a
# priori that a factor has essentially no active gene, which is hard to defend
# for a model whose factors are meant to exist. Prefer bool_var_spec_prob above.
# Section 12e computes, from an existing fit and without refitting, how many
# genes would survive either change.
d_0_power <- 1

tol_abs <- 1e-3
tol_rel <- 1e-5
maxit <- ifelse(bool_dry_run, 1000, 10000)

# Size of the dense grid on which the fitted curves are evaluated. The package
# default is 1000. This is the single largest driver of memory and of the size of
# the saved object: res$list_Y_hat, res$list_Y_low and res$list_Y_upp each hold
# one curve of length n_g for every subject and every gene, that is
#   3 * N * p * n_g * 8 bytes,
# which at N = 150, p = 4471 and n_g = 1000 already amounts to about 16 GB, and
# would reach 80 GB at p = 22354. Over a 50-day window, n_g = 200 gives a
# resolution of a quarter of a day, which is far finer than the data support, and
# divides those objects by five. The grid is also used for the orthonormalisation
# (trapezoidal integration and the SVD of the eigenfunction matrix), so it should
# not be reduced much below that.
n_g <- 200

# Number of cores. More is not better here. bayesSYNC_core calls
# parallel::mclapply about forty times per iteration, and several of those calls
# run over 1:Q, that is five tasks, while the calls inside the q and l loops run
# over 1:N. mclapply forks, so every call duplicates the address space, which
# holds list_cp_C_Y of size (K+2) x p per subject. With roughly two thousand
# forks per iteration, the fork cost dominates once the number of cores exceeds
# the number of tasks that actually have work to do. The run at p = 4471 with 50
# requested cores achieved an effective parallelism of 6, with system time at 37
# per cent of user time, which is the signature of that overhead.
#
# This balance shifts with p. The calls over 1:p give 447 tasks per core at
# p = 22354 on 50 cores, and the work per task in the calls over 1:N grows
# linearly in p, whereas the fork cost grows only with the size of the address
# space. The effective parallelism should therefore improve at large p. Rather
# than guess, run the same short job twice at the target p, once with 16 cores
# and once with 50, and compare the effective parallelism printed after the fit.
#
# On the memory side, fewer forks is strictly safer: copy-on-write does not
# survive the garbage collector, which marks objects and thus triggers the copy.
# If the queue allocates memory per core, requesting 50 cores from SLURM while
# setting n_cpus to 16 gives the larger memory allocation with fewer children.
n_cpus <- ifelse(p < 100, 1, 16) # must match --cpus-per-task in slurm/submit_transcriptomics.sh

# Scaling. Keep this at TRUE. The scaling divides each gene by the standard
# deviation across subjects of its per-subject mean, and the slab on the loadings
# is fixed at N(0, 1), so the data have to be on a comparable scale for that slab
# to mean anything. Turning it off would also make the selection depend on each
# gene's raw variance, which for log-CPM largely reflects expression level and
# transcript length, that is the technical confounder. And it would not buy
# sparsity: in the update Sigma_q_normal_b = 1 / (mu_q_recip_sigsq_eps * rs_tr_qi
# + 1), the +1 is the slab contribution, and once the likelihood term dominates
# it the ratio mu_q_normal_b^2 / Sigma_q_normal_b is nearly scale-invariant.
bool_scale <- T

bool_save <- T

if (bool_save) {
  res_dir <- paste0(out_dir, "/R1_covid_rnaseq_", ifelse(bool_simulated_data, "simulated_", ""),
                    "anneal_", paste0(anneal, collapse = "-"), "_",
                    ifelse(excl_dead, paste0("excl_dead_", dead_definition, "_"), ""),
                    ifelse(bool_excl_subjects,
                           paste0("excl_", paste(vec_excl_subjects, collapse = "-"), "_"), ""),
                    "min_timepts_", min_timepoints, "_",
                    "filter_", gene_filter,
                    ifelse(gene_filter == "none", "", paste0("_q_", gene_quantile)),
                    ifelse(bool_detection_filter, paste0("_det_", detect_min_prop), "_nodet"),
                    ifelse(bool_cov, paste0("_", paste0(covariates, collapse = "_")), ""),
                    "_N_", N, "_p_", p, "_Q_", Q, "_L_", L, "_N_min_", N_min, "_N_max_",
                    N_max, "_K_", ifelse(is.null(K), "auto", K),
                    ifelse(d_0_power == 1, "", paste0("_d0pow_", d_0_power)),
                    ifelse(bool_var_spec_prob, "_omega_jq", ""),
                    "_tol_abs_", tol_abs, "_tol_rel_",
                    tol_rel, "_maxit_", maxit, "_scale_", bool_scale,
                    "_n_cpus_", n_cpus, "_seed_", seed, "/")
  dir.create(res_dir, showWarnings = FALSE)
  sink(paste(res_dir, "out.txt", sep = ""), append = F, split = T, type = "output")
  sink(file(paste(res_dir, "err.txt", sep = ""), open = "wt"), type = "message")
  # the full object is only written on request, so the slim file is what the
  # reload branch below looks for first
  output_file <- file.path(res_dir, "output.RData")
  output_file_slim <- file.path(res_dir, "output_slim.RData")
} else {
  res_dir <- output_file <- output_file_slim <- NULL
}


cat("\n=== ANALYSIS SET ===\n")
cat("Subjects: ", N, "\n", sep = "")
cat("Named subjects excluded: ",
    if (bool_excl_subjects) paste(vec_excl_subjects, collapse = ", ") else "none", "\n", sep = "")
cat("Genes: ", p, "\n", sep = "")
cat("Total observations: ", sum(n_obs), "\n", sep = "")
cat("Observations per subject: min ", N_min, ", median ", median(n_obs),
    ", max ", N_max, "\n", sep = "")
cat("\nDistribution of the number of observations per subject:\n")
print(table(n_obs))
cat("\nSeverity class:\n")
print(table(df_comb_no_HC$severity[match(names(Y), df_comb_no_HC$subject_id)]))


if (!is.null(output_file_slim) && file.exists(output_file_slim)) {
  
  # the slim file holds res_slim, in which the fitted-curve lists are restricted
  # to the plotted variables; the trajectory figures therefore only work for
  # those, which is what vars_plotted records
  load(output_file_slim)
  res <- res_slim
  source("R/fun_application.R")
  source("R/plot_utils.R")
  
} else if (!is.null(output_file) && file.exists(output_file)) {

  # a completed fit with the same settings is reused rather than refitted;
  # delete its folder to refit, for instance after regenerating the data
  # printed with cat, since messages are diverted to err.txt
  cat("Reusing the existing fit in ", output_file, " instead of refitting. ",
      "Delete this folder to refit.\n", sep = "")
  load(output_file)
  source("R/fun_application.R")
  source("R/plot_utils.R")
  
} else {
  
  show_progress <- F
  
  # d_0_power = 1 reproduces the package default, set_hyper(d_0 = p) called
  # internally when list_hyper is NULL
  list_hyper <- set_hyper(d_0 = p^d_0_power)
  
  runtime_bayesSYNC <- system.time(res <- bayesSYNC(time_obs = time_obs,
                                                    Y = Y, L = L, Q = Q,
                                                    K = K,
                                                    n_g = n_g,
                                                    anneal = anneal,
                                                    list_hyper = list_hyper,
                                                    tol_abs = tol_abs,
                                                    tol_rel = tol_rel,
                                                    maxit = maxit,
                                                    seed = seed,
                                                    n_cpus = n_cpus,
                                                    bool_scale = bool_scale,
                                                    bool_var_spec_prob = bool_var_spec_prob,
                                                    show_factor_ppi_progress = show_progress))
  
  res_original <- res
  
  # sign flipping is left out here: the sign convention of the metabolomics
  # analysis was chosen for display and does not transfer to these factors.
  
  if (bool_scale) { # so that the healthy control band is on the scale of the fitted values
    for (vv in all_var) {
      df_comb[, vv] <- (df_comb[, vv] - res$mean_mean_across_subjects[vv]) / res$sd_mean_across_subjects[vv]
    }
  }
  
}

print("=== RUNTIME ===")
print(runtime_bayesSYNC)

# Runtime in hours. system.time() returns five components, and user.self covers
# the master process only: the time spent in the processes forked by mclapply is
# in user.child. The default print method adds them, so they must be added here
# too, otherwise the effective parallelism comes out below one.
user_total <- runtime_bayesSYNC[["user.self"]] +
  ifelse(is.na(runtime_bayesSYNC[["user.child"]]), 0, runtime_bayesSYNC[["user.child"]])
sys_total <- runtime_bayesSYNC[["sys.self"]] +
  ifelse(is.na(runtime_bayesSYNC[["sys.child"]]), 0, runtime_bayesSYNC[["sys.child"]])

cat("\nRuntime: ", format(runtime_bayesSYNC[["elapsed"]] / 3600, digits = 3), " h elapsed, ",
    format(user_total / 3600, digits = 3), " h user (master and children), ",
    format(sys_total / 3600, digits = 3), " h system\n", sep = "")

# Should approach n_cpus. A value far below it means the parallelised part is a
# small share of the work; a large system share points to fork overhead.
cat("Effective parallelism: ", format(user_total / runtime_bayesSYNC[["elapsed"]], digits = 3),
    "x on ", n_cpus, " requested cores, system time at ",
    format(100 * sys_total / max(user_total, 1e-9), digits = 2), "% of user time\n", sep = "")

# The number of variational iterations, not p or N, is usually what separates a
# real-data run from a simulation of comparable size, so it is the first quantity
# to compare across runs. bayesSYNC returns it as i_iter. Note that the full ELBO
# trace is not returned, only its final value in ELBO_iter, so convergence can
# only be inspected through the printed log.
if (!is.null(res$i_iter)) cat("Variational iterations: ", res$i_iter, "\n", sep = "")
cat("Final ELBO: ", format(res$ELBO_iter, digits = 8), "\n", sep = "")


print(paste0("Factor PPIs: ",
             paste0(format(res$factor_ppi, digits = 2), collapse = ", ")))


## ---------------------------------------------------------------------------
## 12. Output
## ---------------------------------------------------------------------------

set.seed(seed)

vec_col_grey <- rep("white", N)
names(vec_col_grey) <- as.character(df_comb$severity[match(names(Y), df_comb$subject_id)])

vec_col_grey[names(vec_col_grey) %in% "A"] <- "grey90"
vec_col_grey[names(vec_col_grey) %in% "B"] <- "grey80"
vec_col_grey[names(vec_col_grey) %in% "C"] <- "grey65"
vec_col_grey[names(vec_col_grey) %in% "D"] <- "grey45"
vec_col_grey[names(vec_col_grey) %in% "E"] <- "grey5"

disp <- Q %/% 5 + 1 # integer division
if (bool_save) {
  pdf(paste0(res_dir, "/scores_1_vs_2.pdf"), height = 3*disp, width = 2.5*ceiling(Q/disp))
}
par(mfrow = c(disp, ceiling(Q/disp)), pty = "s")
null <- sapply(1:Q, function(q) plot(res$list_Zeta_hat[[q]][,1], res$list_Zeta_hat[[q]][,2],
                                     main = paste0("Scores factor q = ", q),
                                     xlab = "Scores FPC 1", ylab = "Scores FPC 2", pch = 20,
                                     col = vec_col_grey))
if (bool_save) {
  dev.off()
  pdf(paste0(res_dir, "/mean_functions.pdf"), width = 8.5, height = 6)
}
par(mfrow = c(2, 3))
sapply(sample(1:p, 6), function(j) plot_mean_functions(j, res$time_g, mu_g = NULL,
                                                       res$list_mu_hat))
if (bool_save) {
  dev.off()
  pdf(paste0(res_dir, "/eigenfunctions.pdf"), width = 3*L, height = 3*Q)
}
cum_pve_thres <- 0.99
par(mfrow = c(Q, L))
null <- sapply(1:Q, function(q)
  plot_factor_eigenfunctions(q, res$time_g, Psi_g = NULL,
                             res$list_list_Phi_hat,
                             res$list_cumulated_pve,
                             cum_pve_thres = cum_pve_thres))
if (bool_save) {
  dev.off()
}

n_random <- 6
vec_random_i <- sample(1:N, n_random)
vec_random_j <- sample(1:p, n_random)
if (bool_save) {
  pdf(paste0(res_dir, "/traj_i_", paste0(vec_random_i, collapse = "-"), "_j_",
             paste0(vec_random_j, collapse = "-"), ".pdf"),
      width = 10, height = 6)
}
par(mfrow = c(2, 3))
for (id in 1:n_random) {
  random_i <- vec_random_i[id]
  random_j <- vec_random_j[id]
  plot_trajectory(random_i, random_j,
                  time_obs, res$time_g, res$Y,
                  res$list_Y_hat, res$list_Y_low, res$list_Y_upp, days_thres = days_thres,
                  df_comb = df_comb)
}
if (bool_save) {
  dev.off()
  pdf(paste0(res_dir, "/h_i_extreme_subjects.pdf"),
      width = 9, height = 2.6*Q)
}

par(mfrow = c(Q, 4))
for (q in 1:Q) {
  for (l in 1:2) {
    plot_factor_trajectory(res$list_Zeta_hat, q, l,
                           res$time_g,
                           res$list_h_hat, res$list_h_low, res$list_h_upp,
                           days_thres = days_thres, df_comb = df_comb)
  }
}
if (bool_save) {
  dev.off()
}


# Top genes per factor. The variable names already carry the HGNC symbol when
# bool_use_symbols is TRUE, and the table below also reports the Ensembl
# identifier and the gene description so that each hit can be checked.
thres <- 0.5
n_top <- 30         # rows printed in the table
n_var_display <- 4  # genes shown in the trajectory figure

for (q in 1:Q) {
  
  if (any(res$ppi[,q] > thres)) {
    
    top_ppi_q <- res$ppi[,q][res$ppi[,q] > thres]
    top_B_hat_q <- res$B_hat[,q][res$ppi[,q] > thres]
    
    o_magnitude <- order(abs(top_B_hat_q), decreasing = TRUE)
    top_genes_q <- names(top_B_hat_q)[o_magnitude]
    
    gene_id_q <- gene_id_of_label[top_genes_q]
    
    df_top_q <- data.frame("label" = top_genes_q,
                           "gene_id" = gene_id_q,
                           "symbol" = df_genes_sel[gene_id_q, "hgnc_symbol"],
                           "ppi" = top_ppi_q[o_magnitude],
                           "loading" = top_B_hat_q[o_magnitude],
                           stringsAsFactors = FALSE)
    
    cat("\n----- FACTOR ", q, ": ", nrow(df_top_q), " genes with PPI > ", thres,
        " (top ", min(n_top, nrow(df_top_q)), " by loading magnitude) -----\n", sep = "")
    print(head(df_top_q, n_top), row.names = FALSE)
    
    n_show <- min(n_var_display, length(top_genes_q))
    
    if (bool_save) {
      pdf(paste0(res_dir, "/Y_i_extreme_subjects_variables_contributing_to_factor_q_", q,
                 "_top_loading_magnitude.pdf"), width = 9, height = 2.2 * n_show)
    }
    plot_factor_variables_trajectory(res$list_Zeta_hat, res$ppi, q, time_obs,
                                     res$time_g, res$Y, res$list_Y_hat,
                                     res$list_Y_low, res$list_Y_upp, thres_ppi = thres,
                                     df_comb = df_comb,
                                     vec_variables = head(top_genes_q, n_show),
                                     days_thres = days_thres, show_ppi = FALSE)
    if (bool_save) {
      dev.off()
    }
    
  }
}

print("omega_hat")
print(res$omega_hat)
print(res$list_cumulated_pve)


## ---------------------------------------------------------------------------
## 12a. Loadings and factor activation
##      The full loading matrix cannot be shown when p is in the thousands, but
##      the sign and magnitude of the leading contributions can. Three views:
##        (i)   a bar chart per factor of the top loadings, signed;
##        (ii)  a heat map of the same genes across all factors, which shows
##              whether a gene contributes to one factor or to several;
##        (iii) the factor posterior inclusion probabilities.
## ---------------------------------------------------------------------------

n_top_loadings <- 20 # genes shown per factor

# genes ranked by loading magnitude within a factor, restricted to those the
# spike-and-slab actually selects
top_loadings_q <- function(res, q, thres, n_top) {
  sel <- res$ppi[, q] > thres
  if (!any(sel)) return(NULL)
  b <- res$B_hat[sel, q]
  b[order(abs(b), decreasing = TRUE)][1:min(n_top, sum(sel))]
}

## (i) signed bar charts, one panel per active factor

if (bool_save) {
  pdf(paste0(res_dir, "/top_loadings_by_factor.pdf"),
      width = 4.2 * Q, height = 0.22 * n_top_loadings + 1.8)
}
par(mfrow = c(1, Q), mar = c(4.5, 8, 3, 1))
for (q in 1:Q) {
  
  b <- top_loadings_q(res, q, thres, n_top_loadings)
  
  if (is.null(b)) {
    plot.new()
    title(main = paste0("Factor ", q, "\n(no gene with PPI > ", thres, ")"))
    next
  }
  
  b <- rev(b) # largest at the top of the panel
  barplot(b, horiz = TRUE, las = 1, cex.names = 0.6,
          col = ifelse(b > 0, "grey30", "grey70"),
          border = NA, xlab = "Estimated loading",
          main = paste0("Factor ", q, " (PPI = ", format(res$factor_ppi[q], digits = 2), ")"))
  abline(v = 0, col = "black")
  
}
if (bool_save) {
  dev.off()
}

## (ii) heat map of the union of the top genes, across all factors

genes_union <- unique(unlist(lapply(1:Q, function(q) names(top_loadings_q(res, q, thres, n_top_loadings)))))

if (length(genes_union) > 0) {
  
  B_sub <- res$B_hat[genes_union, , drop = FALSE]
  
  # loadings are only interpretable where the gene is selected; the others are
  # blanked out rather than shown as small values
  B_sub[res$ppi[genes_union, , drop = FALSE] <= thres] <- NA
  
  # order the genes by the factor they load on most strongly, then by sign
  ord <- order(apply(abs(B_sub), 1, which.max),
               apply(B_sub, 1, function(x) x[which.max(abs(x))]))
  B_sub <- B_sub[ord, , drop = FALSE]
  
  lim <- max(abs(B_sub), na.rm = TRUE)
  breaks <- seq(-lim, lim, length.out = 101)
  cols <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  
  if (bool_save) {
    pdf(paste0(res_dir, "/loadings_heatmap_top_genes.pdf"),
        width = 1.1 * Q + 3, height = 0.16 * nrow(B_sub) + 1.5)
  }
  par(mfrow = c(1, 1), mar = c(4, 8, 3, 4))
  image(x = 1:Q, y = 1:nrow(B_sub), z = t(B_sub), col = cols, breaks = breaks,
        axes = FALSE, xlab = "Factor", ylab = "",
        main = paste0("Loadings of the top ", n_top_loadings, " genes per factor"))
  axis(1, at = 1:Q, labels = 1:Q)
  axis(2, at = 1:nrow(B_sub), labels = rownames(B_sub), las = 1, cex.axis = 0.5, tick = FALSE)
  box()
  if (bool_save) {
    dev.off()
  }
  
  cat("\nGenes shown in the loading heat map: ", nrow(B_sub), "\n", sep = "")
  cat("Genes loading on more than one factor: ",
      sum(rowSums(!is.na(B_sub)) > 1), "\n", sep = "")
  
}

## (iii) factor posterior inclusion probabilities

if (bool_save) {
  pdf(paste0(res_dir, "/factor_ppi.pdf"), width = 1 * Q + 2, height = 3.5)
}
par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3, 1))
bp <- barplot(res$factor_ppi, names.arg = 1:Q, ylim = c(0, 1),
              col = ifelse(res$factor_ppi > 0.5, "grey30", "grey80"), border = NA,
              xlab = "Factor", ylab = "Posterior inclusion probability",
              main = "Factor activation")
abline(h = 0.5, lty = 2)
text(bp, pmin(res$factor_ppi + 0.05, 0.97), format(res$factor_ppi, digits = 2), cex = 0.7)
# number of genes selected by each factor, for context
mtext(paste0("genes with PPI > ", thres, ": ",
             paste0(colSums(res$ppi > thres), collapse = ", ")), side = 1, line = 3.3, cex = 0.7)
if (bool_save) {
  dev.off()
}


## ---------------------------------------------------------------------------
## 12b. Are extreme scores driven by subjects with few observations?
##      Is the small n_i of the patients with extreme scores typical, or related
##      to the scores themselves? With
##      min_timepoints = 1, a subject with a single observation carries little
##      information of its own: the shared mean function and eigenfunctions are
##      estimated from all subjects, but that subject's scores are shrunk
##      towards the prior. Extreme scores should therefore be associated with
##      larger, not smaller, n_i. A negative correlation below would indicate
##      the opposite and would need to be reported.
## ---------------------------------------------------------------------------

n_obs_Y <- sapply(time_obs, length)[names(Y)]

if (bool_save) {
  pdf(paste0(res_dir, "/abs_scores_vs_n_obs.pdf"), height = 3.2*Q, width = 6.5)
}
par(mfrow = c(Q, 2), pty = "s")
for (q in 1:Q) {
  for (l in 1:2) {
    
    abs_scores <- abs(res$list_Zeta_hat[[q]][names(Y), l])
    ct <- cor.test(n_obs_Y, abs_scores, method = "spearman", exact = FALSE)
    
    plot(jitter(n_obs_Y), abs_scores, pch = 20,
         xlab = "Number of observations", ylab = "|Score|",
         main = paste0("Factor q = ", q, ", FPC ", l),
         sub = paste0("Spearman rho = ", format(ct$estimate, digits = 2),
                      ", p = ", format(ct$p.value, digits = 2)))
    
    cat("Factor ", q, ", FPC ", l, ": Spearman rho between n_i and |score| = ",
        format(ct$estimate, digits = 2), " (p = ", format(ct$p.value, digits = 2), ")\n", sep = "")
    
    # number of observations of the four subjects with the most extreme scores,
    # which are the ones displayed in the trajectory figures
    extreme_subj <- names(sort(abs_scores, decreasing = TRUE))[1:4]
    cat("   n_i of the four most extreme subjects: ",
        paste0(n_obs_Y[extreme_subj], collapse = ", "),
        " (median over all subjects: ", median(n_obs_Y), ")\n", sep = "")
    
  }
}
if (bool_save) {
  dev.off()
}


## ---------------------------------------------------------------------------
## 12c. Are the estimated factors capturing residual batch structure?
##      Section 4 shows that the extraction batch correction leaves residual
##      structure on the analysed samples. What matters for the conclusions is
##      whether that structure ends up in the retained factors. If the scores of
##      the retained factors show no association with the batch, the point can be
##      settled in one sentence. If they do, rerun with bool_cov = TRUE and
##      adjust_for = "extractbatch" as a sensitivity analysis.
##
##      The batch is a sample-level variable and 64 of the 165 subjects span
##      several extraction batches, so the batch of the first sample is used as a
##      subject-level summary. The test is therefore approximate.
## ---------------------------------------------------------------------------

df_first <- df_comb_no_HC[order(df_comb_no_HC$days_from_sx_or_swab_imputed), ]
df_first <- df_first[!duplicated(df_first$subject_id), ]
batch_subj <- droplevels(as.factor(df_first$extractbatch[match(names(Y), df_first$subject_id)]))
sev_subj <- droplevels(as.factor(df_first$severity[match(names(Y), df_first$subject_id)]))

cat("\n--- Association between factor scores and extraction batch ---\n")
cat("Reference: association with severity, on the same scores.\n")

for (q in 1:Q) {
  for (l in 1:2) {
    
    sc <- res$list_Zeta_hat[[q]][names(Y), l]
    p_batch <- kruskal.test(sc ~ batch_subj)$p.value
    p_sev <- kruskal.test(sc ~ sev_subj)$p.value
    
    cat(sprintf("Factor %d, FPC %d: batch p = %.3g   severity p = %.3g   (factor PPI = %.2f)\n",
                q, l, p_batch, p_sev, res$factor_ppi[q]))
    
  }
}

if (bool_save) {
  pdf(paste0(res_dir, "/scores_vs_extractbatch.pdf"), height = 3.2*Q, width = 6.5)
}
par(mfrow = c(Q, 2))
for (q in 1:Q) {
  for (l in 1:2) {
    boxplot(res$list_Zeta_hat[[q]][names(Y), l] ~ batch_subj,
            xlab = "Extraction batch", ylab = "Scores",
            main = paste0("Factor q = ", q, ", FPC ", l))
  }
}
if (bool_save) {
  dev.off()
}


## ---------------------------------------------------------------------------
## 12d. Are the factors global technical axes rather than specific signals?
##      Three symptoms point that way when they occur together: a factor selects
##      a large fraction of the genes, its strong loadings nearly all share one
##      sign, and those genes are enriched for long transcripts. The last one is
##      the signature of an RNA-degradation or 3'-coverage bias, since coverage
##      of long transcripts drops first when RNA quality varies.
##
##      Caveat when reading the length results: transcript length is confounded
##      with gene function, long genes being enriched for chromatin and splicing
##      regulators and short genes for translation machinery. A genuine
##      biological axis opposing transcriptional regulation to translation would
##      produce the same length signature. The enrichment is a warning worth
##      following up, not a proof of artefact.
##
##      Note also that filtering on mean expression is itself length-confounded:
##      voom works on log-CPM, which does not correct for transcript length, so
##      longer transcripts accumulate more reads and appear more expressed.
## ---------------------------------------------------------------------------

thres <- 0.5
gene_length <- df_genes_sel[gene_id_of_label[rownames(res$ppi)], "Length"]
names(gene_length) <- rownames(res$ppi)
log_length <- log10(gene_length)

cat("\n--- Sparsity, sign and transcript length of the loadings ---\n")
cat("Median transcript length, analysed genes: ",
    format(median(gene_length, na.rm = TRUE), digits = 5), " bp\n", sep = "")

for (q in 1:Q) {
  
  sel_q <- res$ppi[, q] > thres
  b_q <- res$B_hat[sel_q, q]
  
  # length of the genes with the strongest loadings, against all analysed genes
  top_q <- names(sort(abs(b_q), decreasing = TRUE))[1:min(100, length(b_q))]
  med_top <- median(gene_length[top_q], na.rm = TRUE)
  pct_top <- mean(gene_length < med_top, na.rm = TRUE)
  
  # association over all selected genes, less sensitive to the cut-off
  rho_len <- suppressWarnings(cor(b_q, log_length[sel_q],
                                  method = "spearman", use = "complete.obs"))
  
  cat(sprintf(paste0("Factor %d: %d / %d genes selected (%.0f%%), signs %d positive / %d negative, ",
                     "top-100 median length %.0f bp (percentile %.0f), ",
                     "rho(loading, log length) = %.2f\n"),
              q, sum(sel_q), nrow(res$ppi), 100 * mean(sel_q),
              sum(b_q > 0), sum(b_q < 0), med_top, 100 * pct_top, rho_len))
  
}

# Per-sample length-bias index: within each sample, the correlation between the
# deviation of each gene from its mean and the log transcript length. A sample
# with degraded RNA is depleted in long transcripts and gets a negative value.
# If this index is associated with the factor scores, the factors are picking up
# the technical axis.
E_analysed <- t(as.matrix(df_comb_no_HC[, all_var, drop = FALSE]))
colnames(E_analysed) <- df_comb_no_HC$sample_id_v2

dev_mat <- E_analysed - rowMeans(E_analysed)
lb <- log_length[rownames(dev_mat)] - mean(log_length[rownames(dev_mat)], na.rm = TRUE)
ok <- !is.na(lb)

length_bias <- colSums(dev_mat[ok, ] * lb[ok]) /
  sqrt(colSums(dev_mat[ok, ]^2) * sum(lb[ok]^2))

cat("\nPer-sample length-bias index: range ",
    format(min(length_bias), digits = 2), " to ", format(max(length_bias), digits = 2),
    ", sd ", format(sd(length_bias), digits = 2), "\n", sep = "")

# subject-level summary, to be compared with the subject-level scores
lb_subj <- tapply(length_bias, df_comb_no_HC$subject_id, mean)[names(Y)]

# Reference: the same scores against an inflammation marker, so that a technical
# association can be weighed against a biological one on the same scale.
crp_subj <- NULL
if ("hsCRP_mg_l" %in% names(df_comb_no_HC)) {
  crp_subj <- tapply(df_comb_no_HC$hsCRP_mg_l, df_comb_no_HC$subject_id,
                     function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))[names(Y)]
}

cat("\nSpearman correlation of the factor scores with the length-bias index",
    ifelse(is.null(crp_subj), "", " and with the peak hsCRP"), ":\n", sep = "")
for (q in 1:Q) {
  for (l in 1:2) {
    
    sc <- res$list_Zeta_hat[[q]][names(Y), l]
    ct_lb <- suppressWarnings(cor.test(lb_subj, sc, method = "spearman", exact = FALSE))
    
    if (is.null(crp_subj)) {
      cat(sprintf("Factor %d, FPC %d: length bias rho = %5.2f (p = %.2g)\n",
                  q, l, ct_lb$estimate, ct_lb$p.value))
    } else {
      ct_crp <- suppressWarnings(cor.test(crp_subj, sc, method = "spearman", exact = FALSE))
      cat(sprintf("Factor %d, FPC %d: length bias rho = %5.2f (p = %.2g)   hsCRP rho = %5.2f (p = %.2g)\n",
                  q, l, ct_lb$estimate, ct_lb$p.value, ct_crp$estimate, ct_crp$p.value))
    }
    
  }
}

if (bool_save) {
  pdf(paste0(res_dir, "/scores_vs_length_bias.pdf"), height = 3.2*Q, width = 6.5)
}
par(mfrow = c(Q, 2), pty = "s")
for (q in 1:Q) {
  for (l in 1:2) {
    plot(lb_subj, res$list_Zeta_hat[[q]][names(Y), l], pch = 20,
         xlab = "Length-bias index", ylab = "Scores",
         main = paste0("Factor q = ", q, ", FPC ", l))
  }
}
if (bool_save) {
  dev.off()
}


## ---------------------------------------------------------------------------
## 12e. Would a stronger sparsity prior change the selection?
##      This answers the question without refitting. In the variational update
##      the posterior inclusion probability satisfies
##          logit(PPI_jq) = log BF_jq + Delta_q,
##      where log BF_jq is the evidence carried by gene j for factor q and
##          Delta_q = E[log omega_q] - E[log(1 - omega_q)]
##      is the term contributed by the sparsity prior. Inverting this gives the
##      implied log BF from the fitted PPIs, and the selection can then be
##      recomputed under any other value of Delta_q.
##
##      With omega_q ~ Beta(1, d_0), moving from d_0 = p to d_0 = p^2 shifts
##      Delta_q down by roughly log(p), that is about 8 units here, and pins
##      omega_q near zero, which shifts it further. The grid below spans that
##      range so that the gain can be read off directly.
##
##      Two limits of this calculation. It holds the rest of the variational
##      solution fixed, whereas a genuine refit would let the loadings and the
##      factors readjust, so it is a first-order answer. And PPIs numerically
##      equal to 0 or 1 give infinite logits, so those genes are reported
##      separately as insensitive to the prior.
## ---------------------------------------------------------------------------

delta_grid <- c(0, -2, -4, -8, -12, -16) # shift applied to Delta_q

cat("\n--- Sensitivity of the selection to the sparsity prior ---\n")
cat("Current setting: d_0 = p^", d_0_power, ", bool_var_spec_prob = ", bool_var_spec_prob, ".\n", sep = "")
cat("A shift of about -", format(log(p), digits = 3),
    " corresponds either to d_0 = p^", d_0_power + 1,
    " or, from the shared setting, to bool_var_spec_prob = TRUE, which moves the\n",
    "threshold from roughly log(d_1) - log(c_1) to log(p).\n", sep = "")

cat("\nDistribution of the fitted PPIs, by factor:\n")
for (q in 1:Q) {
  cat(sprintf("Factor %d: ", q))
  print(table(cut(res$ppi[, q], c(-Inf, 0.5, 0.9, 0.99, 1 - 1e-10, Inf),
                  labels = c("<=0.5", "0.5-0.9", "0.9-0.99", "0.99-1", "=1"))))
}

cat("\nNumber of genes with PPI > 0.5 under a shifted prior log-odds:\n")
cat(sprintf("%-10s", "shift"))
cat(sprintf("%10s", paste0("factor ", 1:Q)), "\n", sep = "")

for (d in delta_grid) {
  cat(sprintf("%-10s", format(d)))
  for (q in 1:Q) {
    ppi_q <- pmin(pmax(res$ppi[, q], 1e-12), 1 - 1e-12) # guard against 0 and 1
    logit_ppi <- log(ppi_q / (1 - ppi_q))
    cat(sprintf("%10d", sum(plogis(logit_ppi + d) > 0.5)))
  }
  cat("\n")
}

cat("\nGenes with a fitted PPI numerically equal to 1 (insensitive to the prior):\n")
print(colSums(res$ppi >= 1 - 1e-10))


## ---------------------------------------------------------------------------
## 12f. Is Q too small, or already too large?
##      res$factor_ppi cannot answer this. It is
##          1 - prod_j (1 - ppi_jq),
##      the probability that at least one gene loads on the factor, which
##      saturates at 1 as soon as p is in the thousands: even a uniform ppi of
##      0.001 over 4471 genes gives 1 - exp(-4.47) = 0.99. Every factor showing
##      a PPI of 1 is therefore expected and carries no information about Q.
##
##      What does carry information is redundancy. If an added factor merely
##      splits an existing axis, its scores will correlate with those of another
##      factor and its selected genes will overlap heavily. The two tables below
##      measure both. Strong correlation between the leading scores of two
##      factors, or a Jaccard index close to one between their gene sets, is the
##      signal that Q has gone past what the data support.
## ---------------------------------------------------------------------------

cat("\n--- Redundancy between factors ---\n")

mat_sc <- sapply(1:Q, function(q) res$list_Zeta_hat[[q]][names(Y), 1])
colnames(mat_sc) <- paste0("F", 1:Q)

cat("\nSpearman correlation between the leading scores (FPC 1):\n")
cor_sc <- cor(mat_sc, method = "spearman")
print(round(cor_sc, 2))

cat("\nLargest absolute off-diagonal correlation: ",
    format(max(abs(cor_sc[upper.tri(cor_sc)])), digits = 2), "\n", sep = "")

sets <- lapply(1:Q, function(q) rownames(res$ppi)[res$ppi[, q] > thres])
jac <- outer(1:Q, 1:Q, Vectorize(function(a, b) {
  length(intersect(sets[[a]], sets[[b]])) / length(union(sets[[a]], sets[[b]]))
}))
dimnames(jac) <- list(paste0("F", 1:Q), paste0("F", 1:Q))

cat("\nJaccard index between the selected gene sets:\n")
print(round(jac, 2))

# A gene set that is a large fraction of the analysed genes makes the
# over-representation test powerless, since the target then covers most of the
# background. This is why factors selecting most genes return no enrichment.
cat("\nShare of the analysed genes selected by each factor:\n")
print(round(colMeans(res$ppi > thres), 2))
cat("Over-representation testing has essentially no power above about 0.5.\n")


## ---------------------------------------------------------------------------
## 12g. Shape of the eigenfunctions
##      Whether the leading eigenfunctions carry readable temporal structure,
##      that is whether they are neither flat (a factor constant in time) nor
##      erratic.
##
##      Two points govern how these numbers should be read.
##
##      Only l = 1 and l = 2 are worth judging. Successive eigenfunctions are
##      constrained to be orthogonal to the preceding ones, so their number of
##      sign changes necessarily grows: an oscillating third or fourth component
##      is a property of any functional PCA, not a finding.
##
##      A flat eigenfunction and an erratic one both fail this test, for
##      opposite reasons. The measures below separate the two: the amplitude
##      ratio says whether the curve varies at all, and the number of turning
##      points together with the roughness say whether it varies readably. A
##      leading eigenfunction with one or two turning points over the window has
##      a shape that can be described in words, which is what the argument needs.
## ---------------------------------------------------------------------------

cat("\n--- Shape of the eigenfunctions ---\n")
cat("Only l = 1 and l = 2 are interpretable; higher components must oscillate.\n\n")

eigen_shape <- function(phi, time_g) {
  
  d1 <- diff(phi) / diff(time_g)
  d2 <- diff(d1) / diff(time_g[-1])
  
  # turning points: sign changes of the first derivative
  s <- sign(d1)
  s <- s[s != 0]
  n_turn <- sum(diff(s) != 0)
  
  # roughness: integrated squared second derivative, scaled by the amplitude so
  # that it does not simply track the size of the curve
  amp <- diff(range(phi))
  rough <- sqrt(mean(d2^2)) / ifelse(amp > 0, amp, NA)
  
  # is the curve flat? amplitude relative to its own root mean square
  flat <- amp / sqrt(mean(phi^2))
  
  c("turning_points" = n_turn, "roughness" = rough, "amplitude_ratio" = flat)
  
}

for (q in 1:Q) {
  
  pve_q <- res$list_cumulated_pve[[q]]
  
  for (l in 1:min(2, L)) {
    
    sh <- eigen_shape(res$list_list_Phi_hat[[q]][, l], res$time_g)
    
    cat(sprintf(paste0("Factor %d, eigenfunction %d: cum PVE %5.1f%%, %d turning point(s), ",
                       "roughness %.2f, amplitude ratio %.2f\n"),
                q, l, pve_q[l], sh[["turning_points"]],
                sh[["roughness"]], sh[["amplitude_ratio"]]))
    
  }
  
}

cat("\nA leading eigenfunction with one or two turning points can be described in\n",
    "words; beyond four, the shape is unlikely to carry a reading. An amplitude\n",
    "ratio near zero means the curve is flat.\n", sep = "")

# The regularity of these curves is governed by K much more than by the learnt
# smoothing variance, whose half-Cauchy(A = 1e5) prior is close to uninformative
# and is estimated from a median of two observations per subject. If the leading
# eigenfunctions still oscillate, reduce K rather than expecting the penalty to
# do the work.
if (bool_save) {
  pdf(paste0(res_dir, "/leading_eigenfunctions.pdf"), height = 3 * Q, width = 6.5)
}
par(mfrow = c(Q, 2))
for (q in 1:Q) {
  for (l in 1:min(2, L)) {
    plot(res$time_g * days_thres, res$list_list_Phi_hat[[q]][, l], type = "l",
         xlab = "Days from symptom onset", ylab = "Eigenfunction",
         main = paste0("Factor ", q, ", l = ", l, " (",
                       format(res$list_cumulated_pve[[q]][l], digits = 3), "% cum PVE)"))
    abline(h = 0, lty = 2)
  }
}
if (bool_save) {
  dev.off()
}


## ---------------------------------------------------------------------------
## 13. Scores against clinical variables
## ---------------------------------------------------------------------------

bool_boxplots <- F
if (bool_boxplots) {
  
  method <- "t.test"
  
  scores <- NULL
  for (q in 1:Q) {
    scores <- cbind(scores, res$list_Zeta_hat[[q]][,1:2])
    colnames(scores)[(ncol(scores)-1):ncol(scores)] <- paste0("Factor_", q, "_FPC_", 1:2)
  }
  vec_var <- colnames(scores)
  
  df_scores <- data.frame(scores, "severity" = names(vec_col_grey))
  res_tests <- plot_boxplots(df_scores,
                             vec_var = vec_var,
                             group_var = "severity",
                             vec_var_to_log = NULL,
                             method = method,
                             main = "Scores FPC by severity class",
                             xlab = "Severity class",
                             ylab = "Scores",
                             vec_col = vec_col_grey,
                             save_path = res_dir,
                             bool_log10_disp = FALSE,
                             bool_overall = T,
                             width = 4.6,
                             height = 3.7)
  
  df_scores <- data.frame(scores,
                          "gender" = as.character(df_comb$gender[match(names(Y), df_comb$subject_id)]))
  res_tests <- plot_boxplots(df_scores,
                             vec_var = vec_var,
                             group_var = "gender",
                             vec_var_to_log = NULL,
                             method = method,
                             main = "Scores FPC by gender",
                             xlab = "Gender",
                             ylab = "Scores",
                             vec_col = vec_col_grey,
                             save_path = res_dir,
                             bool_log10_disp = FALSE,
                             bool_overall = F,
                             width = 4.6,
                             height = 3.7)
  
  if (bool_clinical) {
    
    df_scores <- data.frame(scores, "dead" = ifelse(df_status_Y$dead, "dead", "alive"))
    res_tests <- plot_boxplots(df_scores,
                               vec_var = vec_var,
                               group_var = "dead",
                               vec_var_to_log = NULL,
                               method = method,
                               main = "Scores FPC by vital status",
                               xlab = "Vital status",
                               ylab = "Scores",
                               vec_col = vec_col_grey,
                               save_path = res_dir,
                               bool_log10_disp = FALSE,
                               bool_overall = F,
                               width = 4.6,
                               height = 3.7)
  }
  
}


## ---------------------------------------------------------------------------
## 14. Gene set enrichment
##     Same logic as the metabolite version in covid_study.R: GSEA on the ranked
##     loading magnitudes, and an over-representation test on the genes with
##     PPI > 0.5 against the analysed genes as background.
## ---------------------------------------------------------------------------

bool_enrichment <- TRUE
min_pathway_size <- 5
max_pathway_size <- 200 # broad gene sets are dropped: they tend to reach significance
# through their size rather than through biological coherence

# Choice of gene set library. Trade-offs:
#  - "hallmark": 50 MSigDB hallmark sets, built by aggregating and de-duplicating
#    many overlapping sets. Compact, non-redundant and readable, which is what is
#    needed to name a latent factor. Weak point: coarse, it will not resolve a
#    specific pathway. Default here.
#  - "reactome": finely curated, strong coverage of immune signalling, and a good
#    complement to hallmark for the detail behind a hallmark hit. Weak point:
#    many small overlapping sets, so results need manual pruning.
#  - "kegg": consistent with the metabolite analysis of covid_study.R, which
#    helps when comparing the two applications. Weak point: narrow coverage of
#    immune processes and infrequent updates.
#  - "go_bp": exhaustive but highly redundant, typically returning long lists of
#    nested terms. Hard to summarise in a paper.
pathway_db <- c("hallmark", "reactome", "kegg", "go_bp")[1]

if (bool_enrichment) {
  
  require(clusterProfiler)
  require(org.Hs.eg.db)
  if (pathway_db == "hallmark") require(msigdbr)
  if (pathway_db == "reactome") require(ReactomePA)
  
  # Ensembl to Entrez mapping, needed by all four libraries. The model variables
  # are display labels, so the mapping is built from the Ensembl identifiers and
  # a "label" column is carried along to index res$B_hat and res$ppi.
  df_map <- suppressMessages(bitr(selected_genes, fromType = "ENSEMBL",
                                  toType = c("ENTREZID", "SYMBOL"),
                                  OrgDb = org.Hs.eg.db))
  # de-duplicate on both identifiers. Dropping duplicated Ensembl keys is not
  # enough: several Ensembl genes commonly map to the same Entrez gene, and
  # fgsea, which GSEA() calls, rejects a ranking whose names repeat with
  # "Duplicate values in names(stats) not allowed". The error was caught by the
  # tryCatch around run_gsea() below and surfaced only as "No gene set
  # returned.", for every factor of every run, while the over-representation
  # test was unaffected because enricher() takes a set rather than a ranking.
  df_map <- df_map[!duplicated(df_map$ENSEMBL), , drop = FALSE]
  df_map <- df_map[!is.na(df_map$ENTREZID) & !duplicated(df_map$ENTREZID), , drop = FALSE]
  df_map$label <- gene_label_of_id[df_map$ENSEMBL]
  rownames(df_map) <- df_map$label
  
  cat("Genes mapped to Entrez identifiers: ", nrow(df_map), " / ", length(all_var), "\n", sep = "")
  
  background_entrez <- df_map$ENTREZID
  
  # TERM2GENE table for the libraries handled through the generic interface
  if (pathway_db == "hallmark") {
    df_t2g <- msigdbr(species = "Homo sapiens", category = "H")
    df_t2g <- df_t2g[, c("gs_name", "entrez_gene")]
    names(df_t2g) <- c("term", "gene")
    df_t2g$gene <- as.character(df_t2g$gene)
  } else {
    df_t2g <- NULL
  }
  
  # Wrappers so that the two analyses below stay readable
  run_gsea <- function(rank_list) {
    if (pathway_db == "hallmark") {
      GSEA(geneList = rank_list, TERM2GENE = df_t2g, scoreType = "pos",
           minGSSize = min_pathway_size, maxGSSize = max_pathway_size, pvalueCutoff = 1)
    } else if (pathway_db == "reactome") {
      gsePathway(geneList = rank_list, organism = "human", scoreType = "pos",
                 minGSSize = min_pathway_size, maxGSSize = max_pathway_size, pvalueCutoff = 1)
    } else if (pathway_db == "kegg") {
      gseKEGG(geneList = rank_list, organism = "hsa", scoreType = "pos",
              minGSSize = min_pathway_size, maxGSSize = max_pathway_size, pvalueCutoff = 1)
    } else {
      gseGO(geneList = rank_list, OrgDb = org.Hs.eg.db, ont = "BP", scoreType = "pos",
            minGSSize = min_pathway_size, maxGSSize = max_pathway_size, pvalueCutoff = 1)
    }
  }
  
  run_ora <- function(target_entrez) {
    if (pathway_db == "hallmark") {
      enricher(gene = target_entrez, universe = background_entrez, TERM2GENE = df_t2g,
               minGSSize = min_pathway_size, maxGSSize = max_pathway_size, pvalueCutoff = 1)
    } else if (pathway_db == "reactome") {
      enrichPathway(gene = target_entrez, universe = background_entrez, organism = "human",
                    minGSSize = min_pathway_size, maxGSSize = max_pathway_size, pvalueCutoff = 1)
    } else if (pathway_db == "kegg") {
      enrichKEGG(gene = target_entrez, universe = background_entrez, organism = "hsa",
                 minGSSize = min_pathway_size, maxGSSize = max_pathway_size, pvalueCutoff = 1)
    } else {
      enrichGO(gene = target_entrez, universe = background_entrez, OrgDb = org.Hs.eg.db,
               ont = "BP", minGSSize = min_pathway_size, maxGSSize = max_pathway_size,
               pvalueCutoff = 1)
    }
  }
  
  cat("=== ENRICHMENT ANALYSIS FOR GENES CONTRIBUTING TO FACTORS ===\n",
      "=== library: ", pathway_db, ", gene set size between ", min_pathway_size,
      " and ", max_pathway_size, " ===\n", sep = "")
  
  for (factor_id in 1:Q) {
    
    cat(paste0("----- FACTOR: ", factor_id, " -----\n",
               "--- factor PPI: ", format(res$factor_ppi[factor_id], digits = 3), " ---\n",
               "---------------------\n"))
    
    for (useGSEA in c(T, F)) {
      
      res_enrich <- NULL
      
      if (useGSEA) {
        
        cat("*GSEA analysis (ranking based on estimated loadings)*\n")
        
        rank_list <- abs(res$B_hat[rownames(df_map), factor_id])
        names(rank_list) <- df_map$ENTREZID
        rank_list <- rank_list[!is.na(rank_list) & !duplicated(names(rank_list))]
        rank_list <- sort(rank_list, decreasing = TRUE)
        stopifnot(!anyDuplicated(names(rank_list)))
        
        no_results <- (length(unique(rank_list)) == 1)
        
        if (!no_results) {
          set.seed(seed)
          res_enrich <- tryCatch(run_gsea(rank_list), error = function(e) NULL)
        }
        
      } else {
        
        cat("\n*Over-representation test with background (selection based on PPI > 0.5)*\n")
        
        thres <- 0.5
        # rownames(df_map) and rownames(res$ppi) are both display labels
        target_labels <- intersect(rownames(df_map), rownames(res$ppi)[res$ppi[, factor_id] > thres])
        target_entrez <- df_map[target_labels, "ENTREZID"]
        
        cat("Target list (PPI > 0.5): ", length(target_entrez), " genes\n", sep = "")
        print(df_map[target_labels, "SYMBOL"])
        
        no_results <- (length(target_entrez) == 0)
        
        if (!no_results) {
          set.seed(seed)
          res_enrich <- tryCatch(run_ora(target_entrez), error = function(e) NULL)
        }
        
      }
      
      if (no_results) {
        cat("Inactive factor: no enrichment analysis performed.\n")
      } else if (is.null(res_enrich) || nrow(as.data.frame(res_enrich)) == 0) {
        cat("No gene set returned.\n")
      } else {
        df_enrichment <- as.data.frame(res_enrich)
        print(df_enrichment[df_enrichment$pvalue < 0.05,
                            intersect(c("Description", "setSize", "GeneRatio", "NES",
                                        "pvalue", "p.adjust"), names(df_enrichment))])
      }
      
    }
    
  }
  
}


## ---------------------------------------------------------------------------
## 15. Saving
##     save.image() writes the whole workspace, including E_all, W_all and the
##     EList, all of which can be read back from the source file, and including
##     res$list_Y_hat, res$list_Y_low and res$list_Y_upp, which together hold
##     3 * N * p * n_g doubles. Those three lists are the reason the saved object
##     reached 30 GB, and they are only ever read for the handful of subject and
##     gene pairs that appear in the trajectory figures.
##
##     Two files are written. A small one holding everything needed to redo the
##     diagnostics and the tables, which loads in seconds and can be moved off
##     the cluster; and, optionally, the full object for reproducing every figure.
## ---------------------------------------------------------------------------

if (bool_save) {
  
  # variables actually shown in the trajectory figures: the top genes of each
  # factor plus the randomly chosen ones
  vars_plotted <- unique(c(
    unlist(lapply(1:Q, function(q) {
      sel_q <- res$ppi[, q] > thres
      if (!any(sel_q)) return(NULL)
      b_q <- res$B_hat[sel_q, q]
      names(sort(abs(b_q), decreasing = TRUE))[1:min(n_var_display, length(b_q))]
    })),
    all_var[vec_random_j]
  ))
  
  # slim copy of res: the three large lists are restricted to those variables
  res_slim <- res
  for (nm in c("list_Y_hat", "list_Y_low", "list_Y_upp")) {
    res_slim[[nm]] <- lapply(res[[nm]], function(x) x[vars_plotted])
  }
  # res$Y is the standardised copy of Y made by bayesSYNC when bool_scale is
  # TRUE; it is dropped, and the Y saved below is the input on the log-CPM
  # scale, so it must be standardised again before plotting it against the
  # fitted curves (see replot_transcriptomics_from_slim.R)
  res_slim$Y <- NULL
  
  size_full <- format(object.size(res), units = "GB")
  size_slim <- format(object.size(res_slim), units = "MB")
  cat("\nSize of res: ", size_full, " full, ", size_slim, " after restricting the\n",
      "fitted-curve lists to the ", length(vars_plotted), " plotted variables.\n", sep = "")
  
  # small file: enough for sections 12a to 12e, the tables and the enrichment
  save(res_slim, Y, time_obs, df_comb_no_HC, df_status, df_genes_sel,
       gene_id_of_label, gene_label_of_id, all_var, selected_genes,
       n_obs, subject_ids, vars_plotted, runtime_bayesSYNC,
       N, p, Q, L, K, n_g, seed, days_thres, tnew,
       gene_filter, gene_quantile, min_timepoints, excl_dead, dead_definition,
       bool_excl_subjects, vec_excl_subjects,
       bool_cov, bool_scale, bool_var_spec_prob, d_0_power, bool_simulated_data,
       file = file.path(res_dir, "output_slim.RData"))
  
  # full object, only if the trajectory figures need to be redone for other genes
  bool_save_full <- TRUE
  if (bool_save_full) {
    save(res, Y, time_obs, df_comb, df_comb_no_HC,
         file = file.path(res_dir, "output.RData"))
  }
  
}

