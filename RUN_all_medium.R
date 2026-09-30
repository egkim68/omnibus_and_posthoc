# RUN: all five families, medium mode. Nothing to edit: start it with source(), as in README.txt.
FAMILIES <- c("correlation", "twogroup", "anova", "chisq", "posthoc"); MODE <- "medium"; COLLECT <- TRUE
source(file.path(dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")), "run_pipeline.R"))
