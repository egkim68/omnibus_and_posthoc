# EXPLORER 1 of 4: correlation tests on wide data
# Base R only. Writes CSV files only. Run: source("explore_correlation.R") in R or RStudio
# (pasting also works). Finished configurations are checkpointed, so an interrupted run
# resumes when you paste the script again with the same settings.
#
# QUESTION: when many variables can each have missing cells, when do correlation tests
# (Pearson or Spearman) lose validity or power, and which handling (listwise deletion,
# pairwise deletion, mean imputation) protects them? Missing PERCENTAGE per cell is the
# usual summary; this explorer tests whether layout (width relative to rows) and outliers
# matter more.
#
# METHOD: instead of a small fixed grid, every configuration is DRAWN AT RANDOM from wide
# ranges (see sample_parameters). Configurations in the TRAIN set are used to find rules;
# a separate TEST set of fresh configurations checks them (analyze_explorer.R).
# Parameters (p_): rows p_n (30-300), variables p_k (2-40), true exchangeable correlation
# p_effect (0 = no association, so rejections are Type I errors; otherwise power), per-cell
# missing rate p_rate (0-0.30), mechanism p_mechanism (MCAR; MAR: depends on the next
# variable; MNAR: depends on the cell's own value), selection strength p_strength,
# outlier share p_outlier and size p_outlier_size (SD units), marginal shape p_shape
# (normal, heavy-tailed, skewed).
# Descriptors (d_) derived from the parameters help the analysis find thresholds.
# The test is for the correlation of the first two variables. Listwise deletion keeps only
# rows complete on ALL variables; pairwise keeps rows complete on the tested pair; mean
# imputation fills the tested pair and keeps all rows. Analyses with fewer than 6 usable rows
# count as failures (reported as failure_rate), not as non-rejections.
# This is a research starting design, not an established novelty claim.

# "pilot": quick check that the script runs (too small to analyse); "medium": 400 configurations,
# usable results in minutes; "full": the larger run (2,000 configurations,
# 200 replications each). Change the text below, or set the environment variable SIM_MODE
# (run_pipeline.R does this for you).
MODE <- Sys.getenv("SIM_MODE", unset = "medium")
N_TRAIN <- switch(MODE, pilot = 60L, medium = 300L, full = 1500L)    # configurations used to DISCOVER rules
N_TEST <- switch(MODE, pilot = 20L, medium = 100L, full = 500L)      # fresh configurations used to CHECK them
N_REPLICATIONS <- switch(MODE, pilot = 60L, medium = 100L, full = 200L)
# Two seed blocks. The exploration (pilot and medium) uses the exploratory seed. The confirmatory FULL run
# uses a DIFFERENT, declared seed (exploratory seed + 730000), so its designs and random streams are ones
# the exploration never saw.
BASE_SEED <- if (MODE == "full") 2031L + 730000L else 2031L
ALPHA_LEVEL <- 0.05
RESUME <- TRUE
# Change OUTPUT_ROOT or set the environment variable SIM_OUTPUT_ROOT.
OUTPUT_ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "N:/paper-simulation-explorer")
OUTPUT_DIR <- file.path(OUTPUT_ROOT, paste0("correlation_", MODE))
CHECKPOINT_DIR <- file.path(OUTPUT_DIR, "checkpoints")


# ---- shared building blocks ---------------------------------------------------------------
# Bernoulli deletion with the intercept calibrated so the mean probability equals `rate`.
calibrated_missing <- function(eta, rate) {
  n <- length(eta)
  if (rate <= 0) return(rep(FALSE, n))
  intercept <- uniroot(function(a) mean(plogis(a + eta)) - rate,
                       c(-max(eta) - 40, -min(eta) + 40), tol = 1e-8)$root
  runif(n) < plogis(intercept + eta)
}

# Cell-wise contamination: each value is shifted by +/- `size` with probability `prob`.
contaminate <- function(x, prob, size) {
  if (prob <= 0) return(x)
  hit <- runif(length(x)) < prob
  size <- rep_len(size, length(x))
  x[hit] <- x[hit] + size[hit] * sample(c(-1, 1), sum(hit), replace = TRUE)
  x
}

# Draws with population mean 0 and SD 1: normal, heavy-tailed (t with 3 df) or skewed
# (lognormal with sigma 0.8, centred and scaled by its population moments).
draw_shape <- function(n, shape) {
  switch(shape, normal = rnorm(n), heavy = rt(n, 3) / sqrt(3),
         skewed = (rlnorm(n, 0, 0.8) - exp(0.32)) / sqrt((exp(0.64) - 1) * exp(0.64)),
         stop("Unknown shape"))
}

row_out <- function(label, p_value, estimate, n_used)
  data.frame(label = label, p_value = p_value, estimate = estimate, n_used = n_used,
             status = "Success", stringsAsFactors = FALSE)

mean_safe <- function(x) if (length(x)) mean(x) else NA_real_

# The analysis labels a configuration should report (a family can override this after this point).
expected_labels <- function(p) ANALYSIS_LABELS


# ---- configuration sampler and design -----------------------------------------------------
ANALYSIS_LABELS <- c("Listwise pearson", "Listwise spearman", "Pairwise pearson", "Pairwise spearman",
                     "Mean_imputed pearson", "Mean_imputed spearman")
sample_parameters <- function(n) {
  data.frame(
    p_n = round(exp(runif(n, log(30), log(300)))),
    p_k = round(exp(runif(n, log(2), log(40)))),
    p_effect = ifelse(runif(n) < 0.4, 0, round(runif(n, 0.1, 0.5), 2)),
    p_rate = round(runif(n, 0, 0.30), 3),
    p_mechanism = sample(c("MCAR", "MAR", "MNAR"), n, replace = TRUE),
    p_strength = round(runif(n, 0.5, 2.5), 2),
    p_outlier = round(ifelse(runif(n) < 0.3, 0, runif(n, 0.01, 0.10)), 3),
    p_outlier_size = round(runif(n, 4, 12), 1),
    p_shape = sample(c("normal", "heavy", "skewed"), n, replace = TRUE),
    stringsAsFactors = FALSE)
}
make_design <- function(p) NULL
describe <- function(p, design) {
  data.frame(d_expected_listwise_rows = p$p_n * (1 - p$p_rate)^p$p_k,
             d_expected_pairwise_rows = p$p_n * (1 - p$p_rate)^2,
             d_width_ratio = p$p_k / p$p_n)
}

rep_function <- function(p, design) {
  n <- p$p_n; k <- p$p_k
  X <- sqrt(1 - p$p_effect) * matrix(draw_shape(n * k, p$p_shape), n, k) +
    sqrt(p$p_effect) * rnorm(n)
  X <- contaminate(X, p$p_outlier, p$p_outlier_size)
  M <- if (p$p_mechanism == "MCAR") matrix(runif(n * k) < p$p_rate, n, k) else {
    driver <- if (p$p_mechanism == "MAR") X[, c(2:k, 1), drop = FALSE] else X
    matrix(calibrated_missing(p$p_strength * as.vector(driver), p$p_rate), n, k)
  }
  Xm <- X; Xm[M] <- NA
  fill <- function(v) { v[is.na(v)] <- mean(v, na.rm = TRUE); v }
  sets <- list(Listwise = which(rowSums(M) == 0), Pairwise = which(!M[, 1] & !M[, 2]),
               Mean_imputed = seq_len(n))
  rows <- list()
  for (handling in names(sets)) {
    idx <- sets[[handling]]
    x <- if (handling == "Mean_imputed") fill(Xm[, 1]) else Xm[idx, 1]
    y <- if (handling == "Mean_imputed") fill(Xm[, 2]) else Xm[idx, 2]
    for (method in c("pearson", "spearman")) {
      label <- paste(handling, method)
      rows[[label]] <- if (length(idx) < 6 || sd(x) == 0 || sd(y) == 0)
        data.frame(label = label, p_value = NA_real_, estimate = NA_real_,
                   n_used = length(idx), status = "Too few rows", stringsAsFactors = FALSE) else {
        res <- suppressWarnings(cor.test(x, y, method = method, exact = FALSE))
        row_out(label, res$p.value, unname(res$estimate), length(idx))
      }
    }
  }
  do.call(rbind, rows)
}


# ---- running and summarising --------------------------------------------------------------
# One row per configuration and analysis label: rejection rate at alpha with its Monte Carlo
# SE, mean estimate where one exists, mean observations used, and the failure rate.
aggregate_reps <- function(rows, expected = ANALYSIS_LABELS) {
  failed <- sum(rows$label == "FAILED")
  rows <- rows[rows$label != "FAILED", , drop = FALSE]
  pieces <- lapply(split(rows, rows$label), function(d) {
    ok <- d$status == "Success"; s <- d[ok, ]; m <- nrow(s)
    rate <- mean_safe(s$p_value < ALPHA_LEVEL)
    data.frame(label = d$label[1], successful = m,
               failure_rate = 1 - m / N_REPLICATIONS, rejection_rate = rate,
               rejection_mcse = if (m > 0) sqrt(rate * (1 - rate) / m) else NA_real_,
               mean_estimate = mean_safe(s$estimate), mean_n_used = mean_safe(s$n_used),
               stringsAsFactors = FALSE)
  })
  # A label with no successful replication still gets a row (failure rate 1), so that
  # configurations where everything fails are recorded, not lost.
  for (label in setdiff(expected, names(pieces)))
    pieces[[label]] <- data.frame(label = label, successful = 0L, failure_rate = 1,
      rejection_rate = NA_real_, rejection_mcse = NA_real_, mean_estimate = NA_real_,
      mean_n_used = NA_real_, stringsAsFactors = FALSE)
  out <- do.call(rbind, pieces); rownames(out) <- NULL
  out$failed_replications <- failed
  out
}

# Fresh configurations for the CHECK set use a different random stream and seed block.
sample_configs <- function(n, seed, offset) {
  set.seed(seed)
  configs <- sample_parameters(n)
  configs <- cbind(data.frame(config_id = seq_len(n),
                              config_seed = BASE_SEED + offset + 100000000L + seq_len(n)),
                   configs)
  configs
}

run_set <- function(set_name, configs, offset) {
  results <- vector("list", nrow(configs))
  for (i in seq_len(nrow(configs))) {
    p <- configs[i, ]
    checkpoint <- file.path(CHECKPOINT_DIR, sprintf("%s_%05d.rds", set_name, i))
    if (RESUME && file.exists(checkpoint)) res <- readRDS(checkpoint) else {
      design <- make_design(p)
      rows <- lapply(seq_len(N_REPLICATIONS), function(b) {
        set.seed(BASE_SEED + offset + p$config_id * 1000L + b)
        tryCatch(rep_function(p, design), error = function(e) data.frame(label = "FAILED",
          p_value = NA_real_, estimate = NA_real_, n_used = NA_real_,
          status = conditionMessage(e), stringsAsFactors = FALSE))
      })
      res <- aggregate_reps(do.call(rbind, rows), expected_labels(p))
      desc <- describe(p, design)
      res <- cbind(p[rep(1L, nrow(res)), , drop = FALSE],
                   desc[rep(1L, nrow(res)), , drop = FALSE], res)
      rownames(res) <- NULL
      temporary <- paste0(checkpoint, ".tmp")
      saveRDS(res, temporary)
      if (file.exists(checkpoint)) unlink(checkpoint)
      if (!file.rename(temporary, checkpoint)) stop("Checkpoint save failed")
    }
    results[[i]] <- res
    if (i %% 25 == 0 || i == nrow(configs))
      cat(sprintf("[%s] configuration %d/%d\n", set_name, i, nrow(configs)))
  }
  out <- do.call(rbind, results); rownames(out) <- NULL
  out$set <- set_name
  out
}

# Change CODE_VERSION whenever editing functions, or use a new OUTPUT_ROOT.
CODE_VERSION <- "1.0.0"
run_explorer <- function() {
  stopifnot(MODE %in% c("pilot", "medium", "full"), N_TRAIN >= 5, N_TEST >= 5, N_REPLICATIONS >= 10)
  cat(sprintf("Running in %s mode: %d train + %d test configurations, %d replications each.\nResults go to: %s\n",
              toupper(MODE), N_TRAIN, N_TEST, N_REPLICATIONS, OUTPUT_DIR))
  if (MODE == "pilot")
    cat("*** PILOT MODE is only a quick check that the script runs. Its results are too small to analyse.\n",
        "*** For usable results set MODE to \"medium\" (minutes) or \"full\" (the real run) and run again.\n\n")
  dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(CHECKPOINT_DIR))
    stop("Cannot create the output folder ", OUTPUT_DIR, ". Check that the drive exists ",
         "or change OUTPUT_ROOT at the top of the script.", call. = FALSE)
  train <- sample_configs(N_TRAIN, BASE_SEED + 11L, 0L)
  test <- sample_configs(N_TEST, BASE_SEED + 22L, 500000000L)
  config <- list(version = CODE_VERSION, mode = MODE, train = N_TRAIN, test = N_TEST,
                 replications = N_REPLICATIONS, seed = BASE_SEED, alpha = ALPHA_LEVEL,
                 train_configs = train, test_configs = test)
  config_file <- file.path(OUTPUT_DIR, "configuration.rds")
  if (file.exists(config_file) && RESUME) {
    old <- readRDS(config_file)
    if (!identical(old, config)) {
      # Same settings but an older code version: the old results cannot be reused, so the folder is
      # ARCHIVED (renamed, never deleted) and a fresh run starts. Any other change stops the run.
      same_settings <- identical(old[setdiff(names(old), "version")], config[setdiff(names(config), "version")])
      if (!same_settings)
        stop("Configuration changed. Use a new OUTPUT_ROOT or set RESUME <- FALSE.")
      archived <- paste0(OUTPUT_DIR, "_archived_v", as.character(old$version))
      if (dir.exists(archived)) archived <- paste0(archived, "_", format(Sys.time(), "%Y%m%d%H%M%S"))
      cat(sprintf("The results in %s come from code version %s, and this script is version %s.\nArchiving the old folder as %s and starting a fresh run.\n\n",
                  OUTPUT_DIR, as.character(old$version), CODE_VERSION, archived))
      if (!file.rename(OUTPUT_DIR, archived))
        stop("Could not rename ", OUTPUT_DIR, ". Move that folder somewhere else by hand and run again.", call. = FALSE)
      dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)
      saveRDS(config, config_file)
    }
  } else saveRDS(config, config_file)
  writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "session_info.txt"))
  out <- list()
  for (set_name in c("train", "test")) {
    configs <- if (set_name == "train") train else test
    offset <- if (set_name == "train") 0L else 500000000L
    out[[set_name]] <- run_set(set_name, configs, offset)
    write.csv(out[[set_name]], file.path(OUTPUT_DIR,
              sprintf("explore_correlation_%s.csv", set_name)), row.names = FALSE)
  }
  cat("\nFinished. Results:", normalizePath(OUTPUT_DIR), "\n")
  if (MODE == "pilot")
    cat("*** This was a PILOT run (quick check only). Set MODE <- \"medium\" or \"full\" and run again.\n")
  else cat("Analyse them with analyze_explorer.R (it finds the rules and checks them on the test set).\n")
  invisible(out)
}
# The whole workflow runs inside run_explorer(), so a problem stops the run cleanly even if
# you paste line by line. The output folder is created if it is missing.
explorer_results <- run_explorer()
