# EXPLORER 4 of 4: chi-square and Fisher tests when part of the answers are missing
# Base R only. Writes CSV files only. Run: source("explore_chisq.R") in R or RStudio
# (pasting also works). Finished configurations are checkpointed, so an interrupted run
# resumes when you paste the script again with the same settings.
#
# QUESTION: a row variable (group) and a column variable (answer) are cross-tabulated, and
# some answers are missing (nonresponse). Screening suggested that nonresponse depending on
# the group alone, or on the answer alone, cannot create false associations, while nonresponse
# depending on both together does. This explorer maps how strong that interaction has to be, in
# which table shapes, and whether chi-square and Fisher differ.
#
# METHOD: every configuration is DRAWN AT RANDOM from wide ranges (see sample_parameters).
# Configurations in the TRAIN set are used to find rules; a separate TEST set of fresh
# configurations checks them (analyze_explorer.R).
# Parameters (p_): rows p_R (2-6), columns p_C (2-8), sample size p_n (50-1000), row
# concentration p_row_alpha (small = very unequal groups), column concentration p_col_alpha
# (small = sparse, skewed answers), nonresponse rate p_rate (0-0.4), mechanism p_mechanism
# (MCAR; Group: depends on the row only; Outcome: depends on the answer only; Joint: depends on
# row and answer TOGETHER, so it contains an interaction), selection strength p_strength, true
# association p_effect (0 = independence, so rejections are Type I errors; otherwise a random
# tilt of the cell probabilities of that size, so read it as association strength).
# The selection scores per row, per column and per cell are drawn once per configuration.
# Descriptors (d_): number of cells, expected smallest cell count after nonresponse, effective
# sample size, and row imbalance.
# Tests: Pearson chi-square (no continuity correction), Fisher (exact for 2x2, simulated
# p-value with 500 draws otherwise), and the rule "chi-square if every expected count is at
# least 5, else Fisher". The estimate reported is Cramer's V of the observed table.
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
BASE_SEED <- if (MODE == "full") 2034L + 730000L else 2034L
ALPHA_LEVEL <- 0.05
RESUME <- TRUE
# Change OUTPUT_ROOT or set the environment variable SIM_OUTPUT_ROOT.
OUTPUT_ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "N:/paper-simulation-explorer")
OUTPUT_DIR <- file.path(OUTPUT_ROOT, paste0("chisq_", MODE))
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
ANALYSIS_LABELS <- c("Pearson chi-square", "Fisher", "Rule (chi-square if expected >= 5)")
sample_parameters <- function(n) {
  data.frame(
    p_R = sample(2:6, n, replace = TRUE),
    p_C = sample(2:8, n, replace = TRUE),
    p_n = round(exp(runif(n, log(50), log(1000)))),
    p_row_alpha = round(exp(runif(n, log(1), log(20))), 2),
    p_col_alpha = round(exp(runif(n, log(0.3), log(3))), 2),
    p_rate = round(runif(n, 0, 0.4), 3),
    p_mechanism = sample(c("MCAR", "Group", "Outcome", "Joint"), n, replace = TRUE),
    p_strength = round(runif(n, 0.5, 3), 2),
    p_effect = ifelse(runif(n) < 0.6, 0, round(runif(n, 0.1, 0.5), 2)),
    stringsAsFactors = FALSE)
}
# Cell probabilities and selection scores are built once per configuration from its own seed.
make_design <- function(p) {
  set.seed(p$config_seed)
  R <- p$p_R; C <- p$p_C
  rowp <- rgamma(R, p$p_row_alpha); rowp <- rowp / sum(rowp)
  colp <- pmax(rgamma(C, p$p_col_alpha) / 1, 1e-6); colp <- colp / sum(colp)
  cellp <- outer(rowp, colp) * exp(p$p_effect * matrix(rnorm(R * C), R, C))
  list(rowp = rowp, colp = colp, cellp = cellp / sum(cellp),
       eta_row = rnorm(R), eta_col = rnorm(C), eta_cell = matrix(rnorm(R * C), R, C))
}
describe <- function(p, design) {
  data.frame(d_cells = p$p_R * p$p_C, d_n_effective = p$p_n * (1 - p$p_rate),
             d_expected_min_cell = p$p_n * (1 - p$p_rate) * min(design$rowp) * min(design$colp),
             d_row_imbalance = max(design$rowp) / min(design$rowp))
}

rep_function <- function(p, design) {
  R <- p$p_R; C <- p$p_C
  cell <- sample.int(R * C, p$p_n, replace = TRUE, prob = as.vector(design$cellp))
  row <- (cell - 1L) %% R + 1L; col <- (cell - 1L) %/% R + 1L
  gone <- if (p$p_mechanism == "MCAR") runif(p$p_n) < p$p_rate else {
    eta <- switch(p$p_mechanism, Group = design$eta_row[row], Outcome = design$eta_col[col],
                  Joint = design$eta_cell[cbind(row, col)], stop("Unknown mechanism"))
    calibrated_missing(p$p_strength * eta, p$p_rate)
  }
  tab <- table(factor(row[!gone], 1:R), factor(col[!gone], 1:C))
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  if (nrow(tab) < 2 || ncol(tab) < 2) stop("Degenerate table")
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
  fisher <- if (all(dim(tab) == 2)) fisher.test(tab)$p.value else
    fisher.test(tab, simulate.p.value = TRUE, B = 500)$p.value
  rule <- if (min(chi$expected) >= 5) chi$p.value else fisher
  used <- sum(tab)
  v <- sqrt(unname(chi$statistic) / (used * (min(dim(tab)) - 1)))
  rbind(row_out("Pearson chi-square", chi$p.value, v, used),
        row_out("Fisher", fisher, v, used),
        row_out("Rule (chi-square if expected >= 5)", rule, v, used))
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
              sprintf("explore_chisq_%s.csv", set_name)), row.names = FALSE)
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
