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

## ---- paths ----
output_path <- "/scratch/users/kriss1/output/boot_expers"
metadata <- fread(file.path(output_path, "metadata.csv")) %>%
  unique() %>%
  head(25)

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

## ---- beta-boxplots ----
experiment_boxplots(mcombined)

## ---- beta-contours ----
experiment_contours(combined)

## ---- beta-histograms ----
error_histograms(mcombined, c("method + N", "V + D"))

## ---- theta-samples ----
theta <- get_truth_data(metadata, "theta", "i")
combined <- get_samples(metadata, "theta", c("iteration", "variable", "k"))  %>%
  full_join(get_bootstraps(metadata, "theta", "i")) %>%
  left_join(theta)
