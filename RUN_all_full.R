# RUN: all five families, full mode. Nothing to edit: start it with source(), as in README.txt.
FAMILIES <- c("correlation", "twogroup", "anova", "chisq", "posthoc"); MODE <- "full"; COLLECT <- TRUE
source(file.path(dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")), "run_pipeline.R"))
