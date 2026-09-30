# EXPLORER (post hoc) : pairwise procedures after ANOVA-type omnibus tests, with missing data
# Base R only. Writes CSV files only. Run: source("explore_posthoc.R") in R or RStudio (pasting also
# works), or through run_pipeline.R. Finished configurations are checkpointed, so an interrupted run
# resumes when you paste the script again with the same settings.
#
# QUESTION: after an omnibus test, which pairs differ? When missingness leaves some groups tiny and
# the design unbalanced, which post hoc procedures keep the familywise error rate (FWER) and the
# false-pair rate under control, how much power do they keep, and which omnibus + post hoc PIPELINES
# stay coherent? Some failures are fixed by choosing another procedure (precision channel) and some
# by no procedure (bias channel: selection that differs between groups).
#
# PROCEDURES (all base R; each checked against an independent implementation, see below):
#   Tukey-Kramer     studentized range with the pooled variance (Kramer's unbalanced version).
#   Games-Howell     studentized range with each pair's Welch-Satterthwaite df. Its critical values
#                    come from R's qtukey for df >= 2 and from numerical integration for df < 2 (which
#                    happens whenever two groups with 2 observations meet); they are cached on a df
#                    grid per number of groups and interpolated (decisions agree with exact p-values
#                    in 99.98% of 6,000 random cases, 100% for df < 2).
#   Welch-Holm       Welch t tests for every pair with Holm's adjustment.
#   Dunn-Holm        Dunn's test on the pooled ranks (no tie correction: the data are continuous)
#                    with Holm's adjustment.
#   PIPELINES        an omnibus test must be significant before its post hoc test is used:
#                    F > Tukey-Kramer, Welch > Games-Howell, Brown-Forsythe > Games-Howell,
#                    Kruskal-Wallis > Dunn-Holm.
# Validation: Welch-Holm and Dunn-Holm p-values agree with scipy/statsmodels and scikit-posthocs to
# about 1e-15; Tukey-Kramer with statsmodels to the 4 decimals it reports; Games-Howell with pingouin
# to about 1e-5 where R's ptukey is defined.
#
# WHAT IS RECORDED per configuration (a "label" is one measure of one procedure):
#   "FWER <procedure>"           configurations with NO true difference (p_effect = 0): the share of
#                                replications in which ANY pair is declared different.
#   "Partial FWER <procedure>"   configurations with a true difference in some groups: the share of
#                                replications in which any TRULY EQUAL pair is declared different.
#   "Power any-pair <procedure>" the share in which at least one truly different pair is detected.
#   "Pair power <procedure>"     the average share of truly different pairs detected (reported as the
#                                mean estimate; it has no rejection rate).
# The null has two flavours, told apart by p_sd_ratio: 1 means IDENTICAL distributions (the null of
# rank-based procedures such as Dunn), above 1 means EQUAL MEANS with unequal spread (the null of
# mean-based procedures; rank-based procedures answer a different question there).
#
# METHOD: every configuration is DRAWN AT RANDOM from wide ranges (see sample_parameters);
# configurations in the TRAIN set are used to find rules and a separate TEST set of fresh
# configurations checks them (analyze_explorer.R). Parameters (p_): groups p_G (3-40), total size p_N
# (60-600), group-size inequality p_size_cv, SD ratio p_sd_ratio (1 in a quarter of the designs),
# size-spread pairing p_pairing, missing rate p_rate (0-0.4), mechanism p_mechanism (MCAR;
# MAR_group: groups lose data at different rates; MNAR: larger or smaller values deleted more often,
# strength p_strength, sign p_direction), outlier share p_outlier and size p_outlier_size, shape p_shape
# (normal, heavy-tailed, skewed), true shift p_effect (0 = no difference; otherwise the shifted groups
# differ by that many SDs) applied to the share p_alt_share of the groups.
# Descriptors (d_) come from the actual randomly built design (expected group sizes after missingness,
# their spread, the number of pairs). Groups with fewer than 2 observations are dropped; fewer than
# 3 usable groups counts as a failure (failure_rate). alpha = 0.05 for every procedure.
# The Selector rules for the omnibus tests (explore_anova.R) are not repeated here: which pipeline to
# route to is decided after looking at these results.
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
BASE_SEED <- if (MODE == "full") 2035L + 730000L else 2035L
ALPHA_LEVEL <- 0.05
RESUME <- TRUE
# Change OUTPUT_ROOT or set the environment variable SIM_OUTPUT_ROOT.
OUTPUT_ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "N:/paper-simulation-explorer")
OUTPUT_DIR <- file.path(OUTPUT_ROOT, paste0("posthoc_", MODE))
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


PROCEDURES <- c("Tukey-Kramer", "Games-Howell", "Welch-Holm", "Dunn-Holm", "F > Tukey-Kramer",
                "Welch > Games-Howell", "Brown-Forsythe > Games-Howell", "Kruskal-Wallis > Dunn-Holm")
ANALYSIS_LABELS <- c(paste("FWER", PROCEDURES), paste("Partial FWER", PROCEDURES),
                     paste("Power any-pair", PROCEDURES), paste("Pair power", PROCEDURES))
# Null configurations report only FWER labels; configurations with a true difference report the other three.
expected_labels <- function(p) if (p$p_effect == 0) paste("FWER", PROCEDURES) else
  c(paste("Partial FWER", PROCEDURES), paste("Power any-pair", PROCEDURES), paste("Pair power", PROCEDURES))

# ---- configuration sampler and design -----------------------------------------------------
sample_parameters <- function(n) {
  effect <- ifelse(runif(n) < 0.5, 0, round(runif(n, 0.8, 2.0), 2))
  data.frame(
    p_G = round(exp(runif(n, log(3), log(40)))),
    p_N = round(exp(runif(n, log(60), log(600)))),
    p_size_cv = round(runif(n, 0, 1), 2),
    p_sd_ratio = round(ifelse(runif(n) < 0.25, 1, exp(runif(n, log(1.5), log(4)))), 3),
    p_pairing = round(runif(n, -1, 1), 2),
    p_rate = round(runif(n, 0, 0.4), 3),
    p_mechanism = sample(c("MCAR", "MAR_group", "MNAR"), n, replace = TRUE),
    p_strength = round(runif(n, 0.5, 2.5), 2),
    p_direction = sample(c(-1, 1), n, replace = TRUE),
    p_outlier = round(ifelse(runif(n) < 0.3, 0, runif(n, 0.01, 0.10)), 3),
    p_outlier_size = round(runif(n, 4, 12), 1),
    p_shape = sample(c("normal", "heavy", "skewed"), n, replace = TRUE),
    p_effect = effect,
    p_alt_share = ifelse(effect > 0, round(runif(n, 0.1, 0.6), 2), 0),
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
  mu_g <- rep(0, G)
  if (p$p_effect > 0) mu_g[sample.int(G, max(1L, min(G - 1L, as.integer(round(p$p_alt_share * G)))))] <- p$p_effect
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
             d_realised_sd_ratio = max(design$sd_g) / min(design$sd_g), d_total_after = sum(e),
             d_n_pairs = p$p_G * (p$p_G - 1) / 2)
}

# ---- studentized range beyond R's limits ---------------------------------------------------
# R's ptukey needs df >= 2, but a Games-Howell pair of two 2-observation groups has df between 1 and 2.
# For df < 2 the CDF is integrated over the distribution of S = sqrt(chi-square(df) / df).
ptukey_general <- function(q, k, df, lower.tail = TRUE) {
  if (df >= 2) return(ptukey(q, k, df, lower.tail = lower.tail))
  f_s <- function(s) exp(log(2) + (df / 2) * log(df / 2) - lgamma(df / 2) + (df - 1) * log(s) - df * s^2 / 2)
  val <- integrate(function(s) f_s(s) * ptukey(q * s, k, Inf), 0, Inf, rel.tol = 1e-8, subdivisions = 200L)$value
  if (lower.tail) val else 1 - val
}
# Critical values of the studentized range on a df grid, cached per number of groups and interpolated.
QCRIT <- new.env()
qcrit_for <- function(k, df) {
  key <- as.character(k)
  if (is.null(QCRIT[[key]])) {
    grid_df <- c(seq(1, 1.95, length.out = 20), exp(seq(log(2), log(2000), length.out = 60)))
    crit <- vapply(grid_df, function(d) if (d >= 2) qtukey(1 - ALPHA_LEVEL, k, d) else
      uniroot(function(q) ptukey_general(q, k, d, lower.tail = FALSE) - ALPHA_LEVEL, c(0.3, 5000), tol = 1e-7, extendInt = "yes")$root, numeric(1))
    QCRIT[[key]] <- list(x = log(grid_df), y = log(crit))
  }
  g <- QCRIT[[key]]
  exp(approx(g$x, g$y, xout = log(pmin(pmax(df, 1), 2000)), rule = 2)$y)
}
PAIRS <- new.env()
pair_index <- function(k) {
  key <- as.character(k)
  if (is.null(PAIRS[[key]])) PAIRS[[key]] <- which(lower.tri(matrix(0, k, k)), arr.ind = TRUE)
  PAIRS[[key]]
}
# Brown-Forsythe (1974) mean-comparison test: F* with numerator df = k - 1 and a Satterthwaite denominator df.
brown_forsythe_p <- function(y, group) {
  n <- as.vector(table(group)); m <- as.vector(tapply(y, group, mean)); v <- as.vector(tapply(y, group, var))
  N <- sum(n); k <- length(n); grand <- sum(n * m) / N
  c_i <- (1 - n / N) * v
  fstar <- sum(n * (m - grand)^2) / sum(c_i)
  df2 <- 1 / sum((c_i / sum(c_i))^2 / (n - 1))
  pf(fstar, k - 1, df2, lower.tail = FALSE)
}

# ---- the four post hoc procedures: which pairs each declares different ---------------------
pair_decisions <- function(y, gf) {
  n <- as.vector(table(gf)); k <- length(n); N <- length(y)
  m <- as.vector(tapply(y, gf, mean)); v <- as.vector(tapply(y, gf, var))
  idx <- pair_index(k); i <- idx[, 1]; j <- idx[, 2]
  d <- abs(m[i] - m[j])
  mse <- sum((n - 1) * v) / (N - k)
  rej_tk <- d / sqrt(mse / 2 * (1 / n[i] + 1 / n[j])) > qtukey(1 - ALPHA_LEVEL, k, N - k)
  a <- v[i] / n[i]; b <- v[j] / n[j]
  df <- (a + b)^2 / (a^2 / (n[i] - 1) + b^2 / (n[j] - 1))
  rej_gh <- d / sqrt((a + b) / 2) > qcrit_for(k, df)
  rej_wh <- p.adjust(2 * pt(-d / sqrt(a + b), df), "holm") < ALPHA_LEVEL
  rb <- as.vector(tapply(rank(y), gf, mean))
  z <- abs(rb[i] - rb[j]) / sqrt(N * (N + 1) / 12 * (1 / n[i] + 1 / n[j]))
  rej_dh <- p.adjust(2 * pnorm(-z), "holm") < ALPHA_LEVEL
  list(i = i, j = j, rej = list(`Tukey-Kramer` = rej_tk, `Games-Howell` = rej_gh, `Welch-Holm` = rej_wh, `Dunn-Holm` = rej_dh))
}

decision_row <- function(label, rejected, used) row_out(label, if (rejected) 0 else 1, NA_real_, used)
estimate_row <- function(label, value, used)
  data.frame(label = label, p_value = NA_real_, estimate = value, n_used = used, status = "Success", stringsAsFactors = FALSE)

rep_function <- function(p, design) {
  G <- p$p_G
  group <- rep(seq_len(G), times = design$n_g)
  sd_obs <- design$sd_g[group]
  y <- design$mu_g[group] + sd_obs * draw_shape(length(group), p$p_shape)
  y <- contaminate(y, p$p_outlier, p$p_outlier_size * sd_obs)
  gone <- switch(p$p_mechanism,
    MCAR = runif(length(y)) < p$p_rate,
    MAR_group = runif(length(y)) < design$rate_g[group],
    MNAR = calibrated_missing(p$p_direction * p$p_strength * y, p$p_rate),
    stop("Unknown mechanism"))
  y <- y[!gone]; group <- group[!gone]
  keep <- as.integer(names(which(table(group) >= 2)))
  y <- y[group %in% keep]; group <- group[group %in% keep]
  if (length(keep) < 3) stop("Fewer than three usable groups")
  gf <- factor(group, levels = keep)
  if (any(as.vector(tapply(y, gf, var)) <= 0)) stop("A group has zero variance")
  used <- length(y)
  dec <- pair_decisions(y, gf)
  mu <- design$mu_g[keep]
  truly_different <- abs(mu[dec$i] - mu[dec$j]) > 1e-12
  p_f <- oneway.test(y ~ gf, var.equal = TRUE)$p.value
  p_w <- oneway.test(y ~ gf)$p.value
  p_bf <- tryCatch(brown_forsythe_p(y, gf), error = function(e) NA_real_)
  p_kw <- kruskal.test(y ~ gf)$p.value
  gate <- function(rej, p_omnibus) if (isTRUE(p_omnibus < ALPHA_LEVEL)) rej else rep(FALSE, length(rej))
  pipes <- c(dec$rej, list(`F > Tukey-Kramer` = gate(dec$rej[["Tukey-Kramer"]], p_f),
                           `Welch > Games-Howell` = gate(dec$rej[["Games-Howell"]], p_w),
                           `Brown-Forsythe > Games-Howell` = gate(dec$rej[["Games-Howell"]], p_bf),
                           `Kruskal-Wallis > Dunn-Holm` = gate(dec$rej[["Dunn-Holm"]], p_kw)))
  rows <- list()
  for (name in PROCEDURES) {
    rej <- pipes[[name]]
    if (p$p_effect == 0) rows[[length(rows) + 1L]] <- decision_row(paste("FWER", name), any(rej), used)
    else {
      if (any(!truly_different)) rows[[length(rows) + 1L]] <- decision_row(paste("Partial FWER", name), any(rej[!truly_different]), used)
      rows[[length(rows) + 1L]] <- decision_row(paste("Power any-pair", name), any(rej[truly_different]), used)
      # Power refers to the truly different pairs that SURVIVE deletion. If none survives in a replication,
      # pair power is undefined (NaN) and the design's mean pair power is NaN; such designs are excluded from
      # power summaries (3 of 200 exploratory alternative designs). Any-pair power counts that case as no rejection.
      rows[[length(rows) + 1L]] <- estimate_row(paste("Pair power", name), mean(rej[truly_different]), used)
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
              sprintf("explore_posthoc_%s.csv", set_name)), row.names = FALSE)
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
