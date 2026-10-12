#!/usr/bin/env Rscript
# Run from the implementation directory. Read compact summaries, not chains.
source("examples/plotting.R")
source("examples/leaf/example.R")
# First argument: directory containing curves.csv (default: bundled values).
# Second argument: directory for the resulting PDF figures.
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args)) args[1L] else "examples/leaf/cached"
output <- if (length(args) >= 2L) args[2L] else "output/leaf/figures"
render_leaf(read.csv(file.path(input, "curves.csv")), read_leaf_observations(), output)
