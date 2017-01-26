#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Experiments with NMF, under different parameterizations and fitting
## techniques. Here, we consider the parameterizations in the parameterizations
## table and then save the resulting fits in a fits/ directory. These are
## visualized separately.
##
## author: kriss1@stanford.edu

## ---- libraries ----
suppressMessages(library("jsonlite"))
suppressMessages(library("SLURMHelpers"))
base_dir <- "../src/sim/nmf"
source(file.path(base_dir, "nmf_utils.R"))
submit_batch <- FALSE ## so we can recompile the knitr document without submitting cluster jobs

## ---- configuration ----
## create the configuration JSON file
config_path <- file.path(base_dir, "config.json")
batch_dir <- file.path(base_dir, "batch")
fits_dir = file.path(base_dir, "..", "..", "..", "data", "fits")
dir.create(batch_dir)
dir.create(fits_dir)

sim_factors <- list(
  "N" = c(100),
  "P" = c(75, 125),
  "zero_inf_prob" = c(0, 0.2)
)
model_factors <- list(
  "inference" = c("gibbs", "vb"),
  "method" = c(
    file.path(base_dir, "..", "..", "stan", "nmf_gamma_poisson.stan"),
    file.path(base_dir, "..", "..",  "stan", "nmf_gamma_poisson_zero.stan")
  )
)

write_configs(
  sim_factors,
  model_factors,
  n_batches = 4,
  config_path = config_path,
  output_dir = fits_dir
)

## ---- submit-jobs ----
## loop over unique values in the "batch" field of the json file
configs <- fromJSON(config_path, simplifyVector = FALSE)
batches <- sapply(configs, function(x) { x$batch })
batch_opts <- list("mem_alloc" = 6000)

for (i in seq_along(unique(batches))) {
  if (!submit_batch) {
    warning("submit_job == FALSE, not running batch.")
    next
  }

  batch_script <- file.path(batch_dir, paste0("batch-", i, ".sbatch"))
  rscript_file <- file.path(base_dir, "src", "nmf_script.R")
  rscript_cmd <- paste("Rscript", rscript_file, config_path, i)

  create_job(
    batch_script,
    paste0("nmf-", i),
    rscript_cmd,
    batch_opts
  )
  system(paste("sbatch", batch_script))
}
