# RUN: ANOVA, full mode. Nothing to edit: start it with source(), as in README.txt.
FAMILIES <- "anova"; MODE <- "full"; COLLECT <- TRUE
source(file.path(dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")), "run_pipeline.R"))
