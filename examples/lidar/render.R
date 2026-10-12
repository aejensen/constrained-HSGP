#!/usr/bin/env Rscript
# Run from the implementation directory. Read compact summaries, not chains.
source("examples/plotting.R")
source("examples/lidar/example.R")
# First argument: directory containing curves.csv (default: bundled values).
# Second argument: directory for the resulting PDF figures.
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args)) args[1L] else "examples/lidar/cached"
output <- if (length(args) >= 2L) args[2L] else "output/lidar/figures"
render_lidar(read.csv(file.path(input, "curves.csv")), read_lidar_observations("examples/lidar"), output)
