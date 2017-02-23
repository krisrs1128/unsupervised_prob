#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Functions to facilitate posterior predictive checks in the antibiotics study.
##
## author: kriss1@stanford.edu

#' Make histograms comparing simulated and true values
#'
#' @param mx [data.frame] A data.frame of the melted original counts data
#' @param m_sim [data.frame] A data.frame of the melted simulation array.
#'   Includes an extra column for the iteration number.
#' @param n_vis [int] The number of simulated samples to show
#' @return p [ggplot] The histogram plot object
compare_histograms <- function(mx, m_sim, n_vis = 4) {

  ## bind simulated and true data, in order to plot
  iter_max <- max(m_sim$iteration)
  hist_data <- rbind(
    cbind(
      mx %>%
      rename(value = truth),
      iteration = NA,
      type = "truth"
    ),
    cbind(
      m_sim %>%
      filter(iteration %in% round(seq(1, iter_max, length.out = n_vis))) %>%
      rename(value = sim_value),
      type = "simulated"
    )
  )

  ## make histograms
  ggplot(hist_data) +
    geom_histogram(aes(x = asinh(value), fill = type), bins = 100) +
    facet_grid(iteration ~ .) +
    scale_fill_manual(values = c("#377eb8", "#4daf4a")) +
    theme(
      panel.border = element_rect(fill = "transparent", size = 0.5)
    )
}
