# Load the public interface into the environment used by source().
# Paths are resolved relative to this file, so callers may work elsewhere.
.hsgp_source_files <- lapply(sys.frames(), function(frame) frame$ofile)
.hsgp_source_files <- Filter(Negate(is.null), .hsgp_source_files)
.hsgp_loader_file <- if (length(.hsgp_source_files)) tail(.hsgp_source_files, 1L)[[1L]] else NULL
if (is.null(.hsgp_loader_file))
  stop("Load this interface with source('path/to/R/load.R').", call. = FALSE)
.hsgp_root <- normalizePath(file.path(dirname(.hsgp_loader_file), ".."),
                           mustWork = TRUE)
for (.hsgp_module in c("model.R", "fit.R"))
  sys.source(file.path(.hsgp_root, "R", .hsgp_module), envir = environment())
rm(.hsgp_loader_file, .hsgp_source_files, .hsgp_module)
