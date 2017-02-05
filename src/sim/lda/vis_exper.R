#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Script for visualizing output from LDA simulation and evaluation pipeline.
## Three main types of views: Boxplots of proportions estimates, across
## configurations, scatterplots of pairs of proportions estimates, and
## histograms of errors.
##
## author: kriss1@stanford.edu

## ---- libraries-boot-expers ----
library("feather")
library("data.table")
library("plyr")
library("dplyr")
library("tidyr")
library("rstan")
library("ggplot2")
library("ggscaffold")
source("../src/sim/lda/vis_utils.R")
base_dir <- "../src/sim/lda"

## ---- paths ----
output_path <- file.path(base_dir, "..", "..", "..", "data", "fits", "lda-sim")
metadata <- fread(file.path(output_path, "metadata.csv")) %>%
  unique()

## ---- beta-samples ----
beta <- get_truth_data(metadata, "beta") %>%
  rename(variable = v)
combined <- get_samples(metadata, "beta", c("iteration", "k", "variable")) %>%
  full_join(get_bootstraps(metadata, "beta")) %>%
  left_join(beta) %>%
  setDT()

## ---- beta-alignment ----
mcombined <- melt_reshaped_samples(combined)
mcombined <- rbind(
  align_posteriors(mcombined %>% filter(method %in% c("vb", "gibbs"))),
  align_bootstraps(mcombined %>% filter(method == "bootstrap"))
)

combined <- mcombined %>%
  gather(type, value, truth, estimate) %>%
  unite(temp, type, dimension) %>%
  spread(temp, value)

## ---- beta-contours-object ----
unique_V <- unique(mcombined$V)
p <- list()
for (i in seq_along(unique_V)) {
  p[[i]] <- experiment_contours(
    combined %>%
    filter_(sprintf("V == %s", unique_V[i]))
  )
}

## ---- betacontours1 ----
p[[1]]

## ---- betacontours2 ----
p[[2]]

## ---- beta-boxplots-object ----
p <- list()
for (i in seq_along(unique_V)) {
  p[[i]] <- experiment_boxplots(
    mcombined %>%
    filter_(sprintf("V == %s", unique_V[i]))
  )
}

## ---- betaboxplot1 ----
p[[1]]

## ---- betaboxplot2 ----
p[[2]]

## ---- betahistograms ----
error_histograms(mcombined, c("method + V", "D + N"))

## ---- theta-samples ----
theta <- get_truth_data(metadata, "theta", "i")
combined <- get_samples(metadata, "theta", c("iteration", "variable", "k"))  %>%
  full_join(get_bootstraps(metadata, "theta", "i")) %>%
  left_join(theta)
