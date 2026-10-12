#!/usr/bin/env Rscript
# Run from the implementation directory. The optional argument selects output.
source("R/load.R")
source("examples/lidar/example.R")
source("examples/lidar/ep_geometry.R")
# Sampling settings below reproduce the manuscript model. Output CSVs are
# read by render.R, while the fit subdirectory retains draws and diagnostics.
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[1L] else "output/lidar"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

prepared <- prepare_lidar_example("examples/lidar")
proposal <- prepare_derivative_ep(prepared$model)
fit <- sample_hsgp(prepared$model, output_dir = file.path(output, "fit"),
  warmup = 400L, samples = 2000L, seed = 77100L,
  initial = proposal$initial, initial_metric = proposal$inverse_metric,
  fixed_initial_metric = TRUE, max_treedepth = 10L, adapt_delta = .9, stepsize = .03)
# Write numerical summaries first so plotting can be repeated without sampling.
write.csv(lidar_plot_summaries(prepared, fit), file.path(output, "curves.csv"), row.names = FALSE)
write.csv(lidar_predictive_checks(prepared, fit),
          file.path(output, "predictive_checks.csv"), row.names = FALSE)
cat("Saved results in", output, "\nInspect the sampler summary before interpreting them.\n")
