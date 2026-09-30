# PAPER FIGURES: a few composite, multi-panel figures instead of dozens of separate images
# Base R only. Reads one ANOVA explorer run (explore_anova_train.csv, explore_anova_test.csv) and the
# files analyze_explorer.R and evaluate_selectors.R wrote to its "analysis" folder, and writes ONE
# image per figure to <run>/paper_figures (created if missing). Independent of the simulation:
# rerun any time, for example after the full run, without touching the results.
#
# Use in R or RStudio: just source() or paste this file. It looks under OUTPUT_ROOT for the ANOVA runs
# (a "full" run is preferred over "medium", and "medium" over "pilot"). To use one folder, set
# RESULTS_DIR below (or SIM_RESULTS_DIR), or from a terminal:
#   Rscript make_paper_figures.R "N:/.../anova_medium"
#
# THE FIGURES (each panel is lettered; captions belong in the manuscript, not on the image):
#   Figure 1  The precision channel: Type I error of each omnibus test against the smallest expected
#             group after missingness, one panel per test, coloured by missingness mechanism.
#   Figure 2  The two channels: the share of harmful designs (Type I error above HARM_LEVEL) for each
#             test when neither, one, or both channels are active, with 95% Wilson intervals.
#   Figure 3  Held-out AUC of four predictors of a harmful design (missing rate alone, fitted model,
#             tree on all parameters, tree on observable parameters only).
#   Figure 4  What drives failure: share of explained deviance by group of parameters.
#   Figure 5  The selector: A) risk against power on held-out designs; B) harmful share by group-size
#             region for the fixed tests and the selector chosen on the training set.
# Figures 1, 2 and 5B pool the TRAIN and TEST configurations for stability (descriptive views);
# Figures 3, 4 and 5A use only the held-out TEST configurations, as the analysis scripts do.
# Figures 3 and 4 were previously one two-panel figure; they are independent analyses (a held-out
# prediction check and a variance decomposition) that only happened to share row labels, so each is
# now its own figure rather than being paired as though they were one measurement viewed two ways.
#
# Version 1.5: Figure 5 shows the FROZEN rule (FROZEN_RULE below) as the tested rule; a threshold
# re-chosen on the train set appears only as a sensitivity marker. The legends of Figures 1 and 5 no
# longer sit on the data, Figure 3 marks an AUC that could not be estimated, and the exact numbers
# behind Figures 2 and 5B are written next to the images (Figure_2_data.csv, Figure_5B_data.csv).
#
# Each figure is drawn at journal double-column width (FIG_WIDTH_IN, 7.2 inches). Text is ASCII-only so
# the script survives any Windows locale. Change PALETTE, thresholds and sizes in the settings below.

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
FIG_FOLDER <- "paper_figures"         # inside the run's folder (created if missing; files are overwritten)
FIG_FORMAT <- "png"                   # "png", "tiff" or "jpeg"
FIG_WIDTH_IN <- 7.2
FIG_DPI <- 300
HARM_LEVEL <- 0.10                    # Type I error above this makes a design "harmful"
PRECISION_LIMIT <- 4                  # "precision channel" active when the smallest expected group is below this
BIAS_MECHANISMS <- c("MNAR")          # mechanisms that create bias between groups (the "bias channel")
FROZEN_RULE <- "Selector BF<5 else Welch"   # frozen before the confirmatory run (PROTOCOL_v1); do not change
TESTS <- c("F test", "Welch ANOVA", "Alexander-Govern", "Brown-Forsythe", "Kruskal-Wallis")
# Figure 1 shows 4 of the 5 in a 2x2 grid for legibility; Alexander-Govern tracks Welch ANOVA closely
# and its numbers are reported in the text (Section 3.1) and Table 3, so nothing is lost by omitting
# its panel here.
FIG1_TESTS <- c("F test", "Welch ANOVA", "Brown-Forsythe", "Kruskal-Wallis")
MECH_COLORS <- c(MCAR = "#0072B2", MAR_group = "#009E73", MNAR = "#D55E00")
CHANNEL_COLORS <- c(Neither = "#7F7F7F", `Precision only` = "#E69F00", `Bias only` = "#CC79A7", Both = "#000000")

# ---- finding the run (same logic as the other scripts) ------------------------------------------------
describe_folder <- function(d) {
  train <- list.files(d, pattern = "^explore_.*_train\\.csv$"); test <- list.files(d, pattern = "^explore_.*_test\\.csv$")
  list(dir = d, name = basename(d), finished = length(train) == 1 && length(test) == 1)
}
find_anova_run <- function(root) {
  if (!dir.exists(root)) stop("Output folder not found: ", root, ". Check that the drive exists, or set OUTPUT_ROOT.", call. = FALSE)
  dirs <- list.dirs(root, recursive = FALSE)
  dirs <- dirs[grepl("^anova_(pilot|medium|full)$", basename(dirs))]
  infos <- Filter(function(i) i$finished, lapply(dirs, describe_folder))
  if (!length(infos)) stop("No finished ANOVA run under ", root, ". Run explore_anova.R first.", call. = FALSE)
  nm <- vapply(infos, function(i) i$name, character(1))
  i <- .pick_mode(sub("^.*_", "", nm))
  if (is.na(i)) stop("No finished ANOVA ", Sys.getenv("SIM_MODE"), " run under ", root, call. = FALSE)
  infos[[i]]$dir
}

# ---- small helpers -------------------------------------------------------------------------------------
wilson <- function(k, n, z = 1.96) {
  if (n == 0) return(c(NA_real_, NA_real_, NA_real_))
  p <- k / n; den <- 1 + z^2 / n; ctr <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(p, max(0, ctr - half), min(1, ctr + half))
}
panel_letter <- function(letter, line = 0.2, cex = 1.25) mtext(letter, side = 3, adj = 0, line = line, font = 2, cex = cex)
shade <- function(col, a) { m <- col2rgb(col); rgb(m[1, ], m[2, ], m[3, ], alpha = a * 255, maxColorValue = 255) }
# A one-row legend drawn item by item, centred in a strip made with plot.new(), so that the key symbols
# and labels never overlap whatever the device or locale. key(i, x, y, w) draws the i-th key in [x, x + w].
row_legend <- function(y, labels, key, cex = 1, key_w = 0.05, gap = 0.012, space = 0.045) {
  plot.window(c(0, 1), c(0, 1)); op <- par(xpd = NA); on.exit(par(op))
  tw <- strwidth(labels, cex = cex); total <- sum(key_w + gap + tw) + space * (length(labels) - 1)
  if (total > 0.98) { cex <- cex * 0.98 / total; tw <- strwidth(labels, cex = cex); total <- sum(key_w + gap + tw) + space * (length(labels) - 1) }
  x <- (1 - total) / 2
  for (i in seq_along(labels)) { key(i, x, y, key_w); text(x + key_w + gap, y, labels[i], adj = c(0, 0.5), cex = cex); x <- x + key_w + gap + tw[i] + space }
}
nice <- function(x) gsub("-", "-", x, fixed = TRUE)

make_figures <- function(results_dir) {
  if (!FIG_FORMAT %in% c("png", "tiff", "jpeg")) stop("FIG_FORMAT must be png, tiff or jpeg", call. = FALSE)
  train_file <- file.path(results_dir, "explore_anova_train.csv"); test_file <- file.path(results_dir, "explore_anova_test.csv")
  if (!file.exists(train_file) || !file.exists(test_file)) stop("No finished ANOVA results in ", results_dir, call. = FALSE)
  mode_label <- sub("^.*_", "", basename(normalizePath(results_dir)))
  tr <- read.csv(train_file, stringsAsFactors = FALSE); te <- read.csv(test_file, stringsAsFactors = FALSE)
  all_sets <- rbind(tr, te); all_sets <- all_sets[all_sets$successful >= 20, ]
  null_all <- all_sets[all_sets$p_effect == 0, ]
  out_dir <- file.path(results_dir, FIG_FOLDER); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out_dir)) stop("Cannot create ", out_dir, call. = FALSE)
  analysis_dir <- file.path(results_dir, "analysis")
  read_analysis <- function(f) { p <- file.path(analysis_dir, f); if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else NULL }
  validation <- read_analysis("analysis_validation.csv"); drivers <- read_analysis("analysis_drivers.csv")
  sel <- read_analysis("selector_thresholds.csv")
  tests <- TESTS[TESTS %in% null_all$label]
  made <- character(0)

  open_dev <- function(name, height) {
    f <- file.path(out_dir, sprintf("%s.%s", name, FIG_FORMAT))
    a <- list(filename = f, width = FIG_WIDTH_IN, height = height, units = "in", res = FIG_DPI)
    if (FIG_FORMAT == "tiff") a$compression <- "lzw"
    if (FIG_FORMAT == "jpeg") a$quality <- 95
    if (isTRUE(capabilities("cairo"))) a$type <- "cairo"
    do.call(get(FIG_FORMAT), a); f
  }

  # ---- Figure 1: the precision channel (2x2 grid; see FIG1_TESTS above) ------------------------------
  f1_tests <- FIG1_TESTS[FIG1_TESTS %in% null_all$label]
  f1 <- open_dev("Figure_1_precision_channel", 9.0)
  layout(matrix(c(1, 1, 2, 3, 4, 5), 3, 2, byrow = TRUE), heights = c(0.30, 1, 1))
  par(oma = c(4.4, 4.7, 0.4, 0.4))
  par(mar = c(0.2, 0.5, 0.2, 0.5)); plot.new()
  # Two legend rows above the panels, so that no legend covers data: the mechanisms, then the reference lines.
  row_legend(0.72, c("MCAR", "Group-dependent (MAR)", "Outcome-driven (MNAR)"), cex = 1.45,
             key = function(i, x, y, w) segments(x, y, x + w, y, col = MECH_COLORS[i], lwd = 4))
  ref_lty <- c(3, 2, 3); ref_col <- c("gray40", "gray40", "gray60")
  row_legend(0.25, c("Nominal 0.05", sprintf("Harmful above %.2f", HARM_LEVEL), sprintf("Precision limit (%d)", PRECISION_LIMIT)), cex = 1.25,
             key = function(i, x, y, w) segments(x, y, x + w, y, col = ref_col[i], lty = ref_lty[i], lwd = 2))
  par(mar = c(2.0, 1.8, 2.9, 0.6), mgp = c(2.5, 0.9, 0), tcl = -0.35, cex.axis = 1.4, las = 1)
  for (i in seq_len(4)) {
    x <- null_all[null_all$label == f1_tests[i], ]
    plot(x$d_min_expected_n, x$rejection_rate, log = "x", xlim = c(0.3, 90), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "",
         pch = 16, cex = 0.65, col = shade(MECH_COLORS[x$p_mechanism], 0.55))
    rect(0.1, -0.05, 200, HARM_LEVEL, col = shade("#000000", 0.05), border = NA)
    abline(h = 0.05, lty = 3, col = "gray40"); abline(h = HARM_LEVEL, lty = 2, col = "gray40")
    abline(v = PRECISION_LIMIT, lty = 3, col = "gray70")
    for (m in names(MECH_COLORS)) {
      s <- x[x$p_mechanism == m, ]
      if (nrow(s) >= 10) { f <- lowess(log(s$d_min_expected_n), s$rejection_rate, f = 0.6); lines(exp(f$x), f$y, col = MECH_COLORS[[m]], lwd = 3.4) }
    }
    axis(1, at = c(0.5, 1, 2, 5, 10, 20, 50), labels = if (i %in% c(3, 4)) c("0.5", "1", "2", "5", "10", "20", "50") else FALSE, cex.axis = 1.4)
    axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels = if (i %in% c(1, 3)) c("0", "0.25", "0.50", "0.75", "1") else FALSE, cex.axis = 1.4); box()
    title(main = f1_tests[i], cex.main = 1.7, font.main = 1, line = 0.8); panel_letter(LETTERS[i], line = 0.7, cex = 1.7)
  }
  mtext("Smallest expected group size after missingness (log scale)", side = 1, outer = TRUE, line = 2.9, cex = 1.45)
  mtext("Type I error rate", side = 2, outer = TRUE, line = 3.1, cex = 1.45, las = 0)
  dev.off(); made <- c(made, f1)

  # ---- Figure 2: the two channels --------------------------------------------------------------------
  bias <- null_all$p_mechanism %in% BIAS_MECHANISMS; prec <- null_all$d_min_expected_n < PRECISION_LIMIT
  null_all$cell <- ifelse(bias & prec, "Both", ifelse(bias, "Bias only", ifelse(prec, "Precision only", "Neither")))
  cell_names <- names(CHANNEL_COLORS)
  ref <- null_all[null_all$label == tests[1], ]; cell_n <- table(factor(ref$cell, cell_names))
  fig2_rows <- list()
  for (i in seq_along(tests)) for (j in seq_along(cell_names)) {
    x <- null_all[null_all$label == tests[i] & null_all$cell == cell_names[j], ]
    k <- sum(x$rejection_rate > HARM_LEVEL); w <- wilson(k, nrow(x))
    fig2_rows[[length(fig2_rows) + 1L]] <- data.frame(test = tests[i], channels = cell_names[j], n_designs = nrow(x), n_harmful = k,
      share_harmful = round(w[1], 4), wilson_low = round(w[2], 4), wilson_high = round(w[3], 4), stringsAsFactors = FALSE)
  }
  write.csv(do.call(rbind, fig2_rows), file.path(out_dir, "Figure_2_data.csv"), row.names = FALSE)
  f2 <- open_dev("Figure_2_two_channels", 3.9)
  par(mar = c(3.6, 8.6, 1.2, 0.8), mgp = c(2.3, 0.7, 0), tcl = -0.3, cex.axis = 1.0, cex.lab = 1.05, las = 1)
  plot(NA, xlim = c(0, 1), ylim = c(length(tests) + 0.5, 0.5), axes = FALSE, xlab = sprintf("Share of designs with Type I error above %.2f", HARM_LEVEL), ylab = "")
  abline(v = seq(0, 1, 0.25), col = "gray90"); abline(h = seq_along(tests), col = "gray95")
  offsets <- c(-0.27, -0.09, 0.09, 0.27); pchs <- c(21, 22, 23, 24)
  for (i in seq_along(tests)) for (j in seq_along(cell_names)) {
    x <- null_all[null_all$label == tests[i] & null_all$cell == cell_names[j], ]
    w <- wilson(sum(x$rejection_rate > HARM_LEVEL), nrow(x)); y <- i + offsets[j]
    if (!is.na(w[1])) { segments(w[2], y, w[3], y, col = CHANNEL_COLORS[[j]], lwd = 1.6); points(w[1], y, pch = pchs[j], bg = CHANNEL_COLORS[[j]], col = "black", cex = 1.15) }
  }
  axis(1, at = seq(0, 1, 0.25), labels = c("0", "25%", "50%", "75%", "100%"), cex.axis = 1.0); axis(2, at = seq_along(tests), labels = tests, tick = FALSE, cex.axis = 1.05); box()
  legend("bottomright", legend = sprintf("%s (n = %d)", cell_names, as.integer(cell_n)), pch = pchs, pt.bg = CHANNEL_COLORS, col = "black", bty = "o", box.col = "gray80", bg = "white", cex = 0.95, title = "Channels active", pt.cex = 1.3)
  dev.off(); made <- c(made, f2)

  # ---- Figure 3: held-out AUC of four predictors of failure (independent of Figure 4) --------------------
  if (!is.null(validation) && !is.null(drivers)) {
    v <- validation[validation$truth %in% "null" & validation$label %in% tests, ]; v <- v[match(tests[tests %in% v$label], v$label), ]
    views <- c("auc_rate_only", "auc_model", "auc_tree", "auc_tree_observable")
    view_cols <- c("#7F7F7F", "#0072B2", "#D55E00", "#009E73"); view_pch <- c(21, 22, 23, 24); view_off <- c(-0.27, -0.09, 0.09, 0.27)
    f3 <- open_dev("Figure_3_auc", 3.9)
    layout(matrix(1:2, 2, 1), heights = c(0.28, 1))
    par(mar = c(0, 0.2, 0.2, 0.2)); plot.new()
    legend("center", legend = c("Missing rate only", "Fitted model", "Tree (all parameters)", "Tree (observable only)"), pch = view_pch, pt.bg = view_cols, col = "black", ncol = 2, bty = "n", cex = 1.05, pt.cex = 1.4)
    par(mar = c(3.6, 8.6, 0.6, 0.6), mgp = c(2.3, 0.7, 0), tcl = -0.3, cex.axis = 1.05, cex.lab = 1.1, las = 1)
    plot(NA, xlim = c(0.4, 1), ylim = c(nrow(v) + 0.5, 0.5), axes = FALSE, xlab = "AUC on held-out test designs", ylab = "")
    abline(v = 0.5, lty = 3, col = "gray40"); abline(v = seq(0.5, 1, 0.1), col = "gray93"); abline(h = seq_len(nrow(v)), col = "gray95")
    for (i in seq_len(nrow(v))) for (k in seq_along(views)) {
      val <- v[[views[k]]][i]; y <- i + view_off[k]
      if (!is.na(val)) points(val, y, pch = view_pch[k], bg = view_cols[k], col = "black", cex = 1.25)
      else if (views[k] == "auc_tree") text(0.405, y, if (is.na(v$auc_tree_observable[i])) "trees: no rule found" else "tree: no rule found", cex = 0.8, adj = 0, col = "gray30")
      else if (views[k] == "auc_tree_observable" && !is.na(v$auc_tree[i])) text(0.405, y, "observable-only tree: no rule found", cex = 0.8, adj = 0, col = "gray30")
    }
    axis(1, at = seq(0.4, 1, 0.1), cex.axis = 1.05); axis(2, at = seq_len(nrow(v)), labels = v$label, tick = FALSE, cex.axis = 1.1); box()
    dev.off(); made <- c(made, f3)

    # ---- Figure 4: what drives failure, grouped by parameter type (independent of Figure 3) --------------
    groups <- c("Group structure", "Spread and pairing", "Missing rate", "Mechanism", "Data shape, outliers")
    gcol <- c("#56B4E9", "#E69F00", "#7F7F7F", "#D55E00", "#CC79A7")
    group_of <- function(term) ifelse(term %in% c("p_G", "p_N", "p_size_cv"), "Group structure", ifelse(term %in% c("p_sd_ratio", "p_pairing"), "Spread and pairing",
                 ifelse(term == "p_rate", "Missing rate", ifelse(term %in% c("p_mechanism", "p_strength", "p_direction"), "Mechanism", "Data shape, outliers"))))
    dr <- drivers[drivers$truth == "null" & drivers$label %in% v$label, ]; dr$group <- group_of(dr$term)
    mat <- sapply(v$label, function(l) sapply(groups, function(g) sum(dr$importance_pct[dr$label == l & dr$group == g], na.rm = TRUE)))
    f4 <- open_dev("Figure_4_drivers", 3.9)
    layout(matrix(1:2, 2, 1), heights = c(0.28, 1))
    par(mar = c(0, 0.2, 0.2, 0.2)); plot.new()
    legend("center", legend = groups, fill = gcol, border = "white", ncol = 3, bty = "n", cex = 1.05)
    par(mar = c(3.6, 8.6, 0.6, 0.8), mgp = c(2.3, 0.7, 0), tcl = -0.3, cex.axis = 1.05, cex.lab = 1.1, las = 1)
    barplot(mat[, rev(seq_len(ncol(mat))), drop = FALSE], horiz = TRUE, col = gcol, border = "white", axes = FALSE, names.arg = colnames(mat)[rev(seq_len(ncol(mat)))],
            las = 1, cex.names = 1.1, xlim = c(0, max(colSums(mat)) * 1.05), xlab = "")
    title(xlab = "Unique share of null deviance (%)", cex.lab = 1.1, line = 2.3)
    axis(1, cex.axis = 1.05); box()
    dev.off(); made <- c(made, f4)
  } else message("Figures 3 and 4 skipped: run analyze_explorer.R first (it writes analysis_validation.csv and analysis_drivers.csv).")

  # ---- Figure 5: the selector (the FROZEN rule is the tested rule) -------------------------------------
  if (!is.null(sel)) {
    if (!"frozen" %in% names(sel)) sel$frozen <- sel$label == FROZEN_RULE        # older selector tables
    fixed <- sel[sel$kind == "fixed test", ]; sels <- sel[sel$kind == "selector", ]
    sels$thr <- as.integer(sub("^Selector BF<([0-9]+) .*$", "\\1", sels$label)); sels$fb <- sub("^.* else ", "", sels$label)
    f5 <- open_dev("Figure_5_selector", 5.3)
    fb_cols <- c(Welch = "#0072B2", F = "#D55E00")
    # the legend of panel A sits in a strip above the panels, so it can never cover a point or a label
    layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE), widths = c(1.1, 1), heights = c(0.11, 1))
    par(oma = c(0, 0, 0, 0), mar = c(0, 0.5, 0, 0.5)); plot.new()
    # panel A shows only the rule family of the paper (Brown-Forsythe below a threshold, else Welch); the
    # "else F" variants and the threshold re-chosen on the training designs stay in selector_thresholds.csv
    lg_pch <- c(16, 5, 18); lg_col <- c("#444444", fb_cols[["Welch"]], fb_cols[["Welch"]]); lg_cex <- c(1.3, 1.1, 1.7)
    row_legend(0.5, c("Fixed test", "Rule at other thresholds", "Frozen rule (BF<5, else Welch)"), cex = 0.9, key_w = 0.02, gap = 0.008, space = 0.04,
               key = function(i, x, y, w) points(x + w / 2, y, pch = lg_pch[i], col = lg_col[i], cex = lg_cex[i]))
    par(mar = c(4.8, 3.9, 1.4, 0.6), mgp = c(3.1, 0.7, 0), tcl = -0.3, cex.axis = 1.0, cex.lab = 1.0, las = 1)
    sels <- sels[sels$fb == "Welch", ]
    xmax <- max(c(fixed$test_share_harmful, sels$test_share_harmful))
    ylo <- min(c(fixed$test_power, sels$test_power)); yhi <- max(c(fixed$test_power, sels$test_power))
    plot(NA, xlim = c(-0.02, xmax + 0.04), ylim = c(ylo - 0.045, yhi + 0.03), xlab = "", ylab = "Average power (not size-adjusted)")
    title(xlab = sprintf("Share of held-out designs\nwith error above %.2f", HARM_LEVEL), line = 3.0)
    for (fb in "Welch") {
      s <- sels[sels$fb == fb, ]; s <- s[order(s$thr), ]; lines(s$test_share_harmful, s$test_power, col = fb_cols[[fb]], lty = 3)
      mark <- ifelse(s$frozen, 18, 5)
      points(s$test_share_harmful, s$test_power, pch = mark, cex = ifelse(s$frozen, 2.0, 1.0), col = if (fb == "Welch") fb_cols[[fb]] else fb_cols[[fb]])
      # label the frozen rule and the ends of the Welch curve (left of the points), and only the loosest end of
      # the F curve (right of it): its other end sits next to Brown-Forsythe, whose label would collide
      lab <- if (fb == "Welch") s$frozen | s$thr %in% range(s$thr) else s$thr == min(s$thr)
      text(s$test_share_harmful[lab], s$test_power[lab], paste0("<", s$thr[lab]), pos = if (fb == "Welch") 2 else 4, cex = 0.8, col = fb_cols[[fb]], offset = 0.5)
    }
    # fixed-test labels: right of the point, except the two rightmost (left) and Brown-Forsythe (below, clear of the selector curves)
    pos_fixed <- ifelse(fixed$test_share_harmful > 0.3, 2, ifelse(fixed$label == "Brown-Forsythe", 1, 4))
    points(fixed$test_share_harmful, fixed$test_power, pch = 16, cex = 1.4, col = "#444444"); text(fixed$test_share_harmful, fixed$test_power, fixed$label, pos = pos_fixed, cex = 0.85, offset = 0.6)
    panel_letter("A", 0.3)
    frozen_label <- sels$label[sels$frozen][1]
    show <- c("Welch ANOVA", "F test", "Brown-Forsythe", frozen_label); show <- show[show %in% null_all$label]
    null_all$bin <- cut(null_all$d_min_expected_n, c(0, 3, 5, 10, Inf), right = FALSE, labels = c("<3", "3 to <5", "5 to <10", ">=10"))
    share <- sapply(show, function(l) tapply(null_all$rejection_rate[null_all$label == l] > HARM_LEVEL, null_all$bin[null_all$label == l], mean))
    counts <- sapply(show, function(l) table(null_all$bin[null_all$label == l]))
    write.csv(data.frame(bin = rep(rownames(share), ncol(share)), label = rep(colnames(share), each = nrow(share)),
                         n_designs = as.vector(counts), share_harmful = round(as.vector(share), 4)), file.path(out_dir, "Figure_5B_data.csv"), row.names = FALSE)
    par(mar = c(4.8, 3.9, 1.4, 0.6), cex.axis = 1.0, cex.lab = 1.0)
    bar_cols <- c("#0072B2", "#D55E00", "#009E73", "#000000")[seq_along(show)]
    barplot(t(share), beside = TRUE, col = bar_cols, border = "white", ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "",
            names.arg = sub(">=", "\u2265", rownames(share)))
    title(xlab = "Smallest expected group\nafter missingness", line = 3.0)
    axis(2, at = seq(0, 1, 0.25), labels = c("0", "25%", "50%", "75%", "100%"), cex.axis = 1.0); box(); panel_letter("B", 0.3)
    mtext(sprintf("Share of designs with error above %.2f", HARM_LEVEL), side = 2, line = 2.6, cex = 1.0, las = 0)
    legend("topright", legend = sub("^Selector (.*)$", "Frozen rule: \\1", show), fill = bar_cols, border = "white", bty = "n", cex = 0.85)
    dev.off(); made <- c(made, f5)
  } else message("Figure 5 skipped: run evaluate_selectors.R first (it writes selector_thresholds.csv).")

  cat(sprintf("ANOVA (%s): %d figure(s) written to %s\n", mode_label, length(made), normalizePath(out_dir)))
  for (f in made) cat("  ", basename(f), "\n")
  invisible(made)
}

run_figures <- function() {
  target <- if (nzchar(RESULTS_DIR)) RESULTS_DIR else find_anova_run(OUTPUT_ROOT)
  cat("Building paper figures from:", basename(target), "\n")
  make_figures(target)
}
# The whole workflow runs inside run_figures(), so a problem stops the run cleanly even if you paste
# line by line. The figure folder is created if it is missing.
paper_figures <- run_figures()
