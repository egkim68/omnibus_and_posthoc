# EXPLORER 3 of 4: one-way ANOVA, wide (many small groups) versus long (few large groups)
# Base R only. Writes CSV files only. Run: source("explore_anova.R") in R or RStudio
# (pasting also works). Finished configurations are checkpointed, so an interrupted run
# resumes when you paste the script again with the same settings.
#
# QUESTION: how do the F test, Welch ANOVA and Kruskal-Wallis behave when a design has many
# or few groups, unequal group sizes, unequal spreads, missing observations, outliers and
# non-normal data at once? In earlier screening, Welch ANOVA failed badly in wide designs
# and missingness made it worse; this explorer asks how general that is, what rule
# separates safe from unsafe designs, and whether two other heteroscedasticity-robust tests
# (Brown-Forsythe and Alexander-Govern) stay valid where Welch fails.
# Tests: F test, Welch ANOVA, Brown-Forsythe (the original 1974 mean-comparison test, numerator
# df = groups - 1), Alexander-Govern, and Kruskal-Wallis. The Brown-Forsythe and Alexander-Govern
# implementations were checked against statsmodels and SciPy (agreement to about 1e-14).
# The "Selector ..." labels are decision rules: in each replication the observed smallest group
# size decides between Brown-Forsythe and Welch (or F). evaluate_selectors.R compares them with
# using one test throughout.
#
# METHOD: every configuration is DRAWN AT RANDOM from wide ranges (see sample_parameters).
# Configurations in the TRAIN set are used to find rules; a separate TEST set of fresh
# configurations checks them (analyze_explorer.R).
# Parameters (p_): groups p_G (3-40), total size p_N (60-600), group-size inequality
# p_size_cv (0-1), max/min SD ratio p_sd_ratio (1-4), size-spread pairing p_pairing (-1 to 1;
# positive means LARGER groups have LARGER spread; negative means the SMALL groups are the
# noisy ones), missing rate p_rate (0-0.4), mechanism p_mechanism (MCAR; MAR_group: groups
# lose data at different rates, spread p_strength with sign p_direction; MNAR: larger or smaller
# values deleted more often), outlier share p_outlier and size p_outlier_size (SD units),
# shape p_shape (normal, heavy-tailed, skewed), true effect p_effect (0 = all group means
# equal, so rejections are Type I errors; otherwise the SD of group means in within-group SD units).
# Descriptors (d_) come from the actual randomly built design: expected group sizes after
# missingness (minimum, mean, share below 5), the correlation between expected group size and
# group spread, and the realised SD ratio. Groups with fewer than 2 observations are dropped;
# fewer than 2 groups left counts as a failure (failure_rate).
# This is a research starting design, not an established novelty claim.

# "pilot": quick check that the script runs (too small to analyse); "medium": 400 configurations,
# the exploratory phase, usable results in minutes; "full": the confirmatory phase (2,000 configurations,
# 200 replications each, about 1.5 hours per family on a desktop computer). Change the text below, or set the environment variable SIM_MODE
# (run_pipeline.R does this for you).
MODE <- Sys.getenv("SIM_MODE", unset = "medium")
N_TRAIN <- switch(MODE, pilot = 60L, medium = 300L, full = 1500L)    # configurations used to DISCOVER rules
N_TEST <- switch(MODE, pilot = 20L, medium = 100L, full = 500L)      # fresh configurations used to CHECK them
N_REPLICATIONS <- switch(MODE, pilot = 60L, medium = 100L, full = 200L)
# Two seed blocks. The exploration (pilot and medium) uses the exploratory seed. The confirmatory FULL run
# uses a DIFFERENT, declared seed (exploratory seed + 730000), so its designs and random streams are ones
# the exploration never saw.
BASE_SEED <- if (MODE == "full") 2033L + 730000L else 2033L
ALPHA_LEVEL <- 0.05
RESUME <- TRUE
# Change OUTPUT_ROOT or set the environment variable SIM_OUTPUT_ROOT.
OUTPUT_ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "N:/paper-simulation-explorer")
OUTPUT_DIR <- file.path(OUTPUT_ROOT, paste0("anova_", MODE))
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
# SELECTORS: rules that pick a test in every replication from the OBSERVED smallest group size
# (empty groups count as size 0): Brown-Forsythe if it is below the threshold, otherwise Welch
# (or the F test). All thresholds are evaluated on the same replications, so the best one can be
# chosen on the TRAIN set and checked on the TEST set (evaluate_selectors.R).
SELECTOR_THRESHOLDS <- c(3L, 4L, 5L, 6L, 8L, 10L, 12L)
SELECTOR_LABELS <- as.vector(t(outer(SELECTOR_THRESHOLDS, c("Welch", "F"),
                                     function(k, f) sprintf("Selector BF<%d else %s", k, f))))
ANALYSIS_LABELS <- c("F test", "Welch ANOVA", "Brown-Forsythe", "Alexander-Govern", "Kruskal-Wallis",
                     SELECTOR_LABELS)
sample_parameters <- function(n) {
  data.frame(
    p_G = round(exp(runif(n, log(3), log(40)))),
    p_N = round(exp(runif(n, log(60), log(600)))),
    p_size_cv = round(runif(n, 0, 1), 2),
    p_sd_ratio = round(exp(runif(n, log(1), log(4))), 3),
    p_pairing = round(runif(n, -1, 1), 2),
    p_rate = round(runif(n, 0, 0.4), 3),
    p_mechanism = sample(c("MCAR", "MAR_group", "MNAR"), n, replace = TRUE),
    p_strength = round(runif(n, 0.5, 2.5), 2),
    p_direction = sample(c(-1, 1), n, replace = TRUE),
    p_outlier = round(ifelse(runif(n) < 0.3, 0, runif(n, 0.01, 0.10)), 3),
    p_outlier_size = round(runif(n, 4, 12), 1),
    p_shape = sample(c("normal", "heavy", "skewed"), n, replace = TRUE),
    p_effect = ifelse(runif(n) < 0.6, 0, round(runif(n, 0.1, 0.4), 2)),
    stringsAsFactors = FALSE)
}
# The group structure is built once per configuration from its own seed.
make_design <- function(p) {
  set.seed(p$config_seed)
  G <- p$p_G
  w <- exp(p$p_size_cv * rnorm(G))
  # p_N is a TARGET total: planned sizes are rounded and raised to at least 3, so with many groups the
  # realized total can exceed p_N (40 groups need at least 120). The realized sizes are what is analysed.
  n_g <- pmax(3L, as.integer(round(p$p_N * w / sum(w))))
  z <- function(x) (x - mean(x)) / sd(x)
  s <- p$p_pairing * z(rank(n_g + runif(G) * 1e-6)) + sqrt(1 - p$p_pairing^2) * z(rank(runif(G)))
  sd_g <- exp(log(p$p_sd_ratio) * (rank(s) / (G + 1) - 0.5))
  mu_g <- p$p_effect * z(rnorm(G))
  miss_w <- exp(p$p_direction * p$p_strength * rnorm(G))
  # p_sd_ratio is a spread parameter r: the realized largest/smallest SD ratio is r^((G-1)/(G+1)).
  # Under MAR_group, p_rate is a TARGET rate: group rates are rescaled to average p_rate and then capped
  # at 0.9, so when the cap binds the expected missing proportion falls below p_rate.
  rate_g <- if (p$p_mechanism == "MAR_group")
    pmin(0.9, p$p_rate * miss_w / (sum(n_g * miss_w) / sum(n_g))) else rep(p$p_rate, G)
  list(n_g = n_g, sd_g = sd_g, mu_g = mu_g, rate_g = rate_g)
}
describe <- function(p, design) {
  e <- design$n_g * (1 - design$rate_g)
  data.frame(d_min_expected_n = min(e), d_mean_expected_n = mean(e), d_share_small = mean(e < 5),
             d_cor_size_sd = if (sd(e) > 0 && sd(design$sd_g) > 0) cor(e, design$sd_g) else 0,
             d_realised_sd_ratio = max(design$sd_g) / min(design$sd_g),
             d_total_after = sum(e))
}

# Brown-Forsythe (1974) mean-comparison test: F* with numerator df = k - 1 and a Satterthwaite
# denominator df.
brown_forsythe_p <- function(y, group) {
  n <- as.vector(table(group)); m <- as.vector(tapply(y, group, mean)); v <- as.vector(tapply(y, group, var))
  N <- sum(n); k <- length(n); grand <- sum(n * m) / N
  c_i <- (1 - n / N) * v
  fstar <- sum(n * (m - grand)^2) / sum(c_i)
  df2 <- 1 / sum((c_i / sum(c_i))^2 / (n - 1))
  pf(fstar, k - 1, df2, lower.tail = FALSE)
}
# Alexander-Govern (1994) test: each group's precision-weighted t is normalised, and the sum of
# squares is referred to a chi-square distribution with k - 1 df.
alexander_govern_p <- function(y, group) {
  n <- as.vector(table(group)); m <- as.vector(tapply(y, group, mean)); s2 <- as.vector(tapply(y, group, var))
  k <- length(n); w <- n / s2; w <- w / sum(w)
  t <- (m - sum(w * m)) / sqrt(s2 / n); v <- n - 1
  a <- v - 0.5; b <- 48 * a^2; cc <- sqrt(a * log1p(t^2 / v))
  z <- cc + (cc^3 + 3 * cc) / b - (4 * cc^7 + 33 * cc^5 + 240 * cc^3 + 855 * cc) / (10 * b^2 + 8 * b * cc^4 + 1000 * b)
  pchisq(sum(z^2), k - 1, lower.tail = FALSE)
}
rep_function <- function(p, design) {
  group <- rep(seq_len(p$p_G), times = design$n_g)
  sd_obs <- design$sd_g[group]
  y <- design$mu_g[group] + sd_obs * draw_shape(length(group), p$p_shape)
  y <- contaminate(y, p$p_outlier, p$p_outlier_size * sd_obs)
  gone <- switch(p$p_mechanism,
    MCAR = runif(length(y)) < p$p_rate,
    MAR_group = runif(length(y)) < design$rate_g[group],
    MNAR = calibrated_missing(p$p_direction * p$p_strength * y, p$p_rate),
    stop("Unknown mechanism"))
  y <- y[!gone]; group <- group[!gone]
  min_observed <- min(tabulate(group, nbins = p$p_G))   # smallest group as the analyst sees it
  keep_groups <- as.integer(names(which(table(group) >= 2)))
  y <- y[group %in% keep_groups]; group <- factor(group[group %in% keep_groups])
  if (nlevels(group) < 2) stop("Fewer than two usable groups")
  used <- length(y)
  p_f <- oneway.test(y ~ group, var.equal = TRUE)$p.value
  p_w <- oneway.test(y ~ group)$p.value
  p_bf <- tryCatch(brown_forsythe_p(y, group), error = function(e) NA_real_)
  p_ag <- tryCatch(alexander_govern_p(y, group), error = function(e) NA_real_)
  p_kw <- kruskal.test(y ~ group)$p.value
  # A test that cannot be computed (for example a group with zero variance) becomes a failure for
  # that test only, not for the whole replication.
  as_row <- function(label, p) if (is.finite(p)) row_out(label, p, NA_real_, used) else
    data.frame(label = label, p_value = NA_real_, estimate = NA_real_, n_used = used,
               status = "Test could not be computed", stringsAsFactors = FALSE)
  rows <- list(as_row("F test", p_f), as_row("Welch ANOVA", p_w), as_row("Brown-Forsythe", p_bf),
               as_row("Alexander-Govern", p_ag), as_row("Kruskal-Wallis", p_kw))
  for (threshold in SELECTOR_THRESHOLDS) for (fallback in c("Welch", "F")) {
    chosen <- if (min_observed < threshold) p_bf else if (fallback == "Welch") p_w else p_f
    rows[[length(rows) + 1L]] <- as_row(sprintf("Selector BF<%d else %s", threshold, fallback), chosen)
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
CODE_VERSION <- "1.2.0"
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
              sprintf("explore_anova_%s.csv", set_name)), row.names = FALSE)
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
