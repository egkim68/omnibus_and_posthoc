# REBUILD: ANOVA, full mode, from the SAVED results. No simulation is run. Nothing to edit: start it with source().
# Redoes only the reporting (frozen-rule selector table, paper figures, hypothesis checks, upload package).
FAMILIES <- "anova"; MODE <- "full"; COLLECT <- TRUE; REANALYSE <- TRUE
source(file.path(dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")), "run_pipeline.R"))
