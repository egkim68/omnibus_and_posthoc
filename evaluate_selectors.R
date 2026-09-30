# EVALUATE THE ANOVA SELECTORS: does a simple rule beat always using one test?
# Base R only. Reads explore_anova_train.csv and explore_anova_test.csv (from explore_anova.R
# version 1.2.0 or later) and writes a small table and one image to the run's "analysis" folder
# (created if missing). Independent of the explorer and of analyze_explorer.R: rerun any time.
#
# Use in R or RStudio: just source() or paste this file. It looks under OUTPUT_ROOT for the ANOVA
# runs (a "full" run is preferred over "medium", and "medium" over "pilot"). To use one folder,
# set RESULTS_DIR below (or SIM_RESULTS_DIR), or from a terminal:
#   Rscript evaluate_selectors.R "N:/.../anova_medium"
#
# THE RULES. In every replication of every random design, a selector looks at the OBSERVED smallest
# group size after missingness (an empty group counts as 0). If it is below the threshold the
# selector uses Brown-Forsythe; otherwise it uses Welch ANOVA ("else Welch") or the F test
# ("else F"). Thresholds 3, 4, 5, 6, 8, 10 and 12 are all evaluated on the same replications.
#
# THE FROZEN RULE (version 1.5). The confirmatory rule was frozen after the exploration as
# FROZEN_RULE below (Brown-Forsythe when the smallest observed group is below 5, otherwise Welch).
# From 1.5 on, the frozen rule is the MAIN result in every mode: it is marked frozen = TRUE, it is the
# filled marker in the figures, and RUN_INFO.txt reports it. The threshold chosen afresh on the train
# set (below) is kept only as a SENSITIVITY row (chosen_on_train = TRUE): in a confirmatory run it
# shows whether the data would now prefer a different threshold, but it is never the tested rule.
# Up to version 1.4 the re-chosen threshold was reported as the main rule; that was a reporting error.
#
# THE HONEST PART. The threshold is chosen ON THE TRAIN CONFIGURATIONS ONLY and is then evaluated on
# the TEST configurations, which it has never seen. The default rule ("power_under_risk_cap") keeps
# the share of harmful designs at or below RISK_CAP and, among thresholds that do, takes the one with
# the most power; if none does, it takes the lowest-risk threshold. The alternative rule
# ("closest_to_nominal") takes the threshold whose average Type I error is closest to 0.05, which
# tends to drift toward always using Brown-Forsythe. A choice at the edge of the threshold grid is
# flagged: widen SELECTOR_THRESHOLDS in explore_anova.R if that happens. The tables show the full threshold grid on both sets, so you can see whether
# the choice is stable. Power is reported for orientation only: it is NOT size-adjusted, so a test
# that rejects too often looks more powerful than it is.
#
# WHAT TO LOOK FOR. In selector_tradeoff.png, the best position is the lower left of the harmful
# share axis (left) with power as high as possible (up). A good selector sits near the fixed
# tests' best trade-off while avoiding their worst failures.

args <- commandArgs(trailingOnly = TRUE)
RESULTS_DIR <- if (length(args) >= 1) args[1] else Sys.getenv("SIM_RESULTS_DIR", unset = "")

# Which run to use: the mode set by the pipeline (SIM_MODE) if any, otherwise full > medium > pilot.
.pick_mode <- function(modes) {
  m <- Sys.getenv("SIM_MODE", unset = "")
  pref <- if (m %in% c("pilot", "medium", "full")) m else c("full", "medium", "pilot")
  i <- which.min(match(modes, pref))
  if (length(i)) i else NA_integer_
}
OUTPUT_ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "N:/paper-simulation-explorer")
ANALYSIS_FOLDER <- "analysis"         # inside the ANOVA run's folder (created if missing)
HARM_LEVEL <- 0.10                    # a null configuration is "harmful" above this Type I error
NOMINAL <- 0.05
SELECTION_RULE <- "power_under_risk_cap"   # or "closest_to_nominal"
RISK_CAP <- 0.10                      # highest acceptable share of harmful designs (train set)
FROZEN_RULE <- "Selector BF<5 else Welch"   # frozen before the confirmatory run (PROTOCOL_v1); do not change
FIXED_TESTS <- c("F test", "Welch ANOVA", "Brown-Forsythe", "Alexander-Govern", "Kruskal-Wallis")
FIG_FORMAT <- "png"                   # "png", "tiff" or "jpeg"
FIG_DPI <- 300

# ---- finding finished runs (same logic as analyze_explorer.R) ------------------------------------
describe_folder <- function(d) {
  train <- list.files(d, pattern = "^explore_.*_train\\.csv$")
  test <- list.files(d, pattern = "^explore_.*_test\\.csv$")
  list(dir = d, name = basename(d), finished = length(train) == 1 && length(test) == 1,
       n_files = length(train) + length(test), checkpoints = length(list.files(file.path(d, "checkpoints"))))
}
explain_folder <- function(info) {
  what <- if (info$finished) "finished" else if (info$checkpoints > 0)
    sprintf("UNFINISHED: %d checkpoints but no results CSV yet; paste the explorer script again with the same settings to resume", info$checkpoints)
  else if (info$n_files > 0) "incomplete results files"
  else "no results yet (wrong folder, or the explorer has not run)"
  sprintf("  %s: %s", info$name, what)
}
find_anova_run <- function(root) {
  if (!dir.exists(root))
    stop("Output folder not found: ", root, ". Check that the drive exists, or set OUTPUT_ROOT ",
         "(or RESULTS_DIR) at the top of the script.", call. = FALSE)
  infos <- Filter(function(i) grepl("^anova_(pilot|medium|full)$", i$name), lapply(list.dirs(root, recursive = FALSE), describe_folder))
  if (!length(infos)) stop("No ANOVA explorer folder under ", root, ". Run explore_anova.R first.", call. = FALSE)
  done <- Filter(function(i) i$finished, infos)
  if (!length(done))
    stop("No finished ANOVA run under ", root, ". What is there:\n",
         paste(vapply(infos, explain_folder, character(1)), collapse = "\n"), call. = FALSE)
  nm <- vapply(done, function(i) i$name, character(1))
  i <- .pick_mode(sub("^.*_", "", nm))
  if (is.na(i)) stop("No finished ANOVA ", Sys.getenv("SIM_MODE"), " run under ", root, call. = FALSE)
  done[[i]]$dir
}

evaluate_selectors <- function(results_dir) {
  if (!FIG_FORMAT %in% c("png", "tiff", "jpeg")) stop("FIG_FORMAT must be png, tiff or jpeg", call. = FALSE)
  if (!dir.exists(results_dir)) stop("Folder not found: ", results_dir, call. = FALSE)
  train_file <- list.files(results_dir, pattern = "^explore_anova_train\\.csv$", full.names = TRUE)
  test_file <- list.files(results_dir, pattern = "^explore_anova_test\\.csv$", full.names = TRUE)
  if (length(train_file) != 1 || length(test_file) != 1)
    stop("No finished ANOVA results in ", results_dir, ".\n", explain_folder(describe_folder(results_dir)), call. = FALSE)
  mode_label <- sub("^.*_", "", basename(normalizePath(results_dir)))
  tr <- read.csv(train_file, stringsAsFactors = FALSE); te <- read.csv(test_file, stringsAsFactors = FALSE)
  selector_labels <- sort(unique(grep("^Selector BF<[0-9]+ else (Welch|F)$", tr$label, value = TRUE)))
  if (!length(selector_labels))
    stop("This run has no 'Selector ...' results. Run explore_anova.R version 1.2.0 or later (in a ",
         "new folder, or after moving the old one away).", call. = FALSE)
  out_dir <- file.path(results_dir, ANALYSIS_FOLDER)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)        # created if missing
  if (!dir.exists(out_dir)) stop("Cannot create ", out_dir, call. = FALSE)

  # One row per configuration, one column per analysis label (its rejection rate).
  wide <- function(d) {
    d <- d[d$successful >= 20, c("config_id", "p_effect", "label", "rejection_rate")]
    w <- reshape(d, idvar = c("config_id", "p_effect"), timevar = "label", direction = "wide")
    names(w) <- sub("^rejection_rate\\.", "", names(w)); w
  }
  w_tr <- wide(tr); w_te <- wide(te)
  metrics <- function(w, label) {
    n <- w[w$p_effect == 0, label]; a <- w[w$p_effect != 0, label]
    c(typeI_mean = mean(n, na.rm = TRUE), share_harmful = mean(n > HARM_LEVEL, na.rm = TRUE), n_harmful = sum(n > HARM_LEVEL, na.rm = TRUE),
      typeI_gap = mean(abs(n - NOMINAL), na.rm = TRUE), power = mean(a, na.rm = TRUE),
      n_null = sum(!is.na(n)), n_alt = sum(!is.na(a)))
  }
  all_labels <- c(FIXED_TESTS[FIXED_TESTS %in% names(w_tr)], selector_labels)
  grid <- do.call(rbind, lapply(all_labels, function(l) {
    a <- metrics(w_tr, l); b <- metrics(w_te, l)
    data.frame(label = l, kind = if (l %in% selector_labels) "selector" else "fixed test",
      train_typeI_mean = a[["typeI_mean"]], train_share_harmful = a[["share_harmful"]], train_typeI_gap = a[["typeI_gap"]], train_power = a[["power"]],
      test_typeI_mean = b[["typeI_mean"]], test_share_harmful = b[["share_harmful"]], test_typeI_gap = b[["typeI_gap"]], test_power = b[["power"]],
      n_null_test = b[["n_null"]], n_alt_test = b[["n_alt"]],
      n_harmful_train = a[["n_harmful"]], n_null_train = a[["n_null"]], n_harmful_test = b[["n_harmful"]],
      all_share_harmful = (a[["n_harmful"]] + b[["n_harmful"]]) / (a[["n_null"]] + b[["n_null"]]),   # train + test pooled ("all designs")
      stringsAsFactors = FALSE)
  }))
  # The threshold for each fallback is chosen on the TRAIN set only.
  chosen <- character(0); at_edge <- character(0)
  for (fallback in c("Welch", "F")) {
    cand <- grid[grepl(paste0(" else ", fallback, "$"), grid$label) & grid$kind == "selector", ]
    thr <- as.integer(sub("^Selector BF<([0-9]+) .*$", "\\1", cand$label))
    pick <- if (SELECTION_RULE == "closest_to_nominal") which.min(cand$train_typeI_gap) else {
      ok <- which(cand$train_share_harmful <= RISK_CAP)
      if (length(ok)) ok[which.max(cand$train_power[ok])] else which.min(cand$train_share_harmful)
    }
    chosen <- c(chosen, cand$label[pick])
    if (thr[pick] %in% range(thr)) at_edge <- c(at_edge, cand$label[pick])
  }
  grid$chosen_on_train <- grid$label %in% chosen
  grid$frozen <- grid$label == FROZEN_RULE
  if (!any(grid$frozen)) stop("The frozen rule '", FROZEN_RULE, "' is not in this run's results.", call. = FALSE)
  grid$role <- ifelse(grid$kind == "fixed test", "fixed test", ifelse(grid$frozen, "FROZEN RULE (tested)",
                 ifelse(grid$chosen_on_train, "re-chosen on train (sensitivity only)", "other threshold")))
  evaluation <- grid[grid$kind == "fixed test" | grid$frozen | grid$chosen_on_train, ]
  evaluation <- evaluation[order(match(evaluation$role, c("FROZEN RULE (tested)", "fixed test", "re-chosen on train (sensitivity only)")), evaluation$test_share_harmful), ]
  num <- vapply(grid, is.numeric, logical(1)); grid[num] <- lapply(grid[num], round, 3)
  num <- vapply(evaluation, is.numeric, logical(1)); evaluation[num] <- lapply(evaluation[num], round, 3)
  write.csv(grid, file.path(out_dir, "selector_thresholds.csv"), row.names = FALSE)
  write.csv(evaluation, file.path(out_dir, "selector_evaluation.csv"), row.names = FALSE)

  # ---- image: harmful share against power on the TEST configurations -----------------------------
  args <- list(filename = file.path(out_dir, sprintf("selector_tradeoff.%s", FIG_FORMAT)),
               width = 10, height = 7.5, units = "in", res = FIG_DPI)
  if (FIG_FORMAT == "tiff") args$compression <- "lzw"
  if (FIG_FORMAT == "jpeg") args$quality <- 95
  if (isTRUE(capabilities("cairo"))) args$type <- "cairo"
  do.call(get(FIG_FORMAT), args)
  par(mar = c(5, 5, 4, 1))
  fixed <- grid[grid$kind == "fixed test", ]; sel <- grid[grid$kind == "selector", ]
  xr <- range(c(fixed$test_share_harmful, sel$test_share_harmful), na.rm = TRUE)
  yr <- range(c(fixed$test_power, sel$test_power), na.rm = TRUE)
  plot(NA, xlim = xr + c(-0.03, 0.06), ylim = yr + c(-0.03, 0.04), xlab = sprintf("Share of test configurations with Type I error above %.2f (lower is better)", HARM_LEVEL),
       ylab = "Average power (not size-adjusted)", main = sprintf("Fixed tests and selectors on fresh test configurations: ANOVA (%s)", mode_label))
  cols <- c(Welch = "#0072B2", F = "#D55E00")
  for (fallback in c("Welch", "F")) {
    s <- sel[grepl(paste0(" else ", fallback, "$"), sel$label), ]
    s <- s[order(as.integer(sub("^Selector BF<([0-9]+) .*$", "\\1", s$label))), ]
    lines(s$test_share_harmful, s$test_power, col = cols[[fallback]], lty = 3)
    points(s$test_share_harmful, s$test_power, pch = ifelse(s$frozen, 18, ifelse(s$chosen_on_train, 9, 5)), cex = ifelse(s$frozen, 2.2, 1.3), col = cols[[fallback]])
    text(s$test_share_harmful, s$test_power, sub("^Selector BF<([0-9]+) .*$", "<\\1", s$label), pos = 3, cex = 0.65, col = cols[[fallback]])
  }
  points(fixed$test_share_harmful, fixed$test_power, pch = 16, cex = 1.6, col = "#444444")
  text(fixed$test_share_harmful, fixed$test_power, fixed$label, pos = 4, cex = 0.8)
  legend("bottomright", c("Fixed test", "Selector, else Welch", "Selector, else F test", "Frozen rule (tested)", "Re-chosen on train (sensitivity)"),
         pch = c(16, 5, 5, 18, 9), col = c("#444444", cols[["Welch"]], cols[["F"]], "black", "black"), bty = "n", cex = 0.85)
  mtext(sprintf("Filled diamond = the frozen rule (%s); crossed diamonds = thresholds re-chosen on the TRAIN set; labels show the threshold", sub("^Selector ", "", FROZEN_RULE)), side = 3, line = 0.2, cex = 0.7)
  dev.off()

  fr <- grid[grid$frozen, ]
  cat(sprintf("\nANOVA selectors (%s). FROZEN RULE %s: held-out harmful share %.3f, all designs %.3f (%d of %d), power %.3f\n",
              mode_label, FROZEN_RULE, fr$test_share_harmful, fr$all_share_harmful, fr$n_harmful_train + fr$n_harmful_test, fr$n_null_train + fr$n_null_test, fr$test_power))
  cat(sprintf("Sensitivity only: rule '%s' would re-choose on train = %s\n", SELECTION_RULE, paste(chosen, collapse = "; ")))
  if (length(at_edge)) cat("NOTE: chosen at the edge of the threshold grid (", paste(at_edge, collapse = "; "), "); widen SELECTOR_THRESHOLDS.\n", sep = "")
  print(evaluation[, c("label", "role", "test_typeI_mean", "test_share_harmful", "all_share_harmful", "test_power")], row.names = FALSE)
  cat("Wrote selector_evaluation.csv, selector_thresholds.csv and selector_tradeoff.", FIG_FORMAT, " to ", normalizePath(out_dir), "\n", sep = "")
  invisible(evaluation)
}

run_evaluation <- function() {
  target <- if (nzchar(RESULTS_DIR)) RESULTS_DIR else find_anova_run(OUTPUT_ROOT)
  cat("Evaluating selectors in:", basename(target), "\n")
  evaluate_selectors(target)
}
# The whole workflow runs inside run_evaluation(), so a problem stops the run cleanly even if you
# paste line by line. The output folder is created if it is missing.
selector_results <- run_evaluation()
