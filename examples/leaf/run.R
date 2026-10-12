#!/usr/bin/env Rscript
# Run from the implementation directory. The optional argument selects output.
source("R/load.R")
source("examples/leaf/example.R")
# Sampling settings below reproduce the manuscript model. Output CSVs are
# read by render.R, while the fit subdirectory retains draws and diagnostics.
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[1L] else "output/leaf"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

prepared <- prepare_leaf_example()
fit <- sample_hsgp(prepared$model, output_dir = file.path(output, "fit"),
  warmup = 1500L, samples = 3000L, seed = 20261010L,
  initial = leaf_initial(prepared$model),
  max_treedepth = 12L, adapt_delta = .99, stepsize = .01)
ordinary <- draw_ordinary_hsgp(prepared$model, samples = 20000L, seed = 20264000L)
# Write numerical summaries first so plotting can be repeated without sampling.
write.csv(leaf_plot_summaries(prepared, fit, ordinary), file.path(output, "curves.csv"), row.names = FALSE)
write.csv(leaf_scalar_summaries(prepared, fit, ordinary), file.path(output, "application_summary.csv"), row.names = FALSE)
write.csv(leaf_predictive_checks(prepared, fit, ordinary), file.path(output, "predictive_checks.csv"), row.names = FALSE)
cat("Saved results in", output, "\nInspect the sampler summary before interpreting them.\n")
