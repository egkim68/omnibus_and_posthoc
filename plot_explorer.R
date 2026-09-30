# PLOTS FOR THE EXPLORERS: one folder per family, one subfolder per test
# Base R only. Reads explore_<family>_train.csv and explore_<family>_test.csv and writes IMAGE
# FILES (PNG by default). Independent of the explorers and of analyze_explorer.R: rerun any time.
#
# Use in R or RStudio: just source() or paste this file. It looks under OUTPUT_ROOT for every
# FINISHED explorer run (a "full" run is preferred over a "pilot" run of the same family) and
# plots each. To plot a single folder instead, set RESULTS_DIR below (or the environment variable
# SIM_RESULTS_DIR), or from a terminal:  Rscript plot_explorer.R "N:/.../anova_full"
#
# WHERE THE IMAGES GO (created if missing; each run OVERWRITES files of the same name, so to keep
# a version, change PLOT_FOLDER first, for example to "plots_v2"):
#   <OUTPUT_ROOT>/<family>_<mode>/plots/<test name>/   one folder per test or handling
#   <OUTPUT_ROOT>/<family>_<mode>/plots/_overview/     comparisons across the tests of a family
#
# PLOTS PER TEST (each point is one random configuration; TRAIN and TEST sets are pooled, see
# USE_SETS; horizontal lines mark the nominal 0.05 and the "harmful" level 0.10):
#   1_typeI_vs_missing_rate     Type I error against the missing percentage, coloured by mechanism.
#   2_typeI_vs_<key design 1>   the same against the first key design variable (see KEY_VARIABLES).
#   3_typeI_vs_<key design 2>   and against the second.
#   4_typeI_map_<1>_by_<2>      average Type I error in bins of the two key variables (a map of
#                               where the test is safe and where it fails).
#   5_power_vs_<key design 1>   power on configurations with a true effect (if there are enough).
#   6_failure_vs_<key design 1> share of replications where the analysis could not run (only if
#                               failures occurred).
# If the missing percentage were a good guide, plot 1 would show a clear trend; the point of plots 2
# to 4 is that a design variable usually explains far more.

args <- commandArgs(trailingOnly = TRUE)
# Leave RESULTS_DIR empty to plot every finished run found under OUTPUT_ROOT.
RESULTS_DIR <- if (length(args) >= 1) args[1] else Sys.getenv("SIM_RESULTS_DIR", unset = "")

# Which run to use: the mode set by the pipeline (SIM_MODE) if any, otherwise full > medium > pilot.
.pick_mode <- function(modes) {
  m <- Sys.getenv("SIM_MODE", unset = "")
  pref <- if (m %in% c("pilot", "medium", "full")) m else c("full", "medium", "pilot")
  i <- which.min(match(modes, pref))
  if (length(i)) i else NA_integer_
}
OUTPUT_ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "N:/paper-simulation-explorer")
PLOT_FOLDER <- "plots"
USE_SETS <- c("train", "test")        # which configuration sets to plot
FIG_FORMAT <- "png"                   # "png", "tiff" or "jpeg"
FIG_WIDTH_IN <- 8
FIG_HEIGHT_IN <- 6
FIG_DPI <- 300
HARM_LEVEL <- 0.10
EXCLUDE_LABELS <- "^(Selector|Partial FWER|Power any-pair|Pair power)"   # rules and power-side post hoc measures are handled elsewhere
# The design variables that matter most for each family (edit freely; p_ = parameter,
# d_ = descriptor derived from the parameters; see the header of each explorer script).
KEY_VARIABLES <- list(
  correlation = c("d_expected_listwise_rows", "p_k"),
  twogroup = c("d_pairing", "p_n_total"),
  anova = c("d_min_expected_n", "p_G"),
  posthoc = c("d_min_expected_n", "p_G"),
  chisq = c("d_expected_min_cell", "p_n"))

# ---- finding finished runs -------------------------------------------------------------------
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
find_runs <- function(root) {
  if (!dir.exists(root))
    stop("Output folder not found: ", root, ". Check that the drive exists, or set OUTPUT_ROOT ",
         "(or RESULTS_DIR) at the top of the script.", call. = FALSE)
  dirs <- list.dirs(root, recursive = FALSE)
  dirs <- dirs[grepl("_(pilot|medium|full)$", basename(dirs)) & !grepl("^_ALL_", basename(dirs))]   # only explorer run folders
  infos <- lapply(dirs, describe_folder)
  if (!length(infos)) stop("No explorer folders found under ", root, ". Run an explorer script first.", call. = FALSE)
  done <- Filter(function(i) i$finished, infos)
  if (!length(done))
    stop("No finished explorer run under ", root, ". What is there:\n",
         paste(vapply(infos, explain_folder, character(1)), collapse = "\n"), call. = FALSE)
  family <- sub("_(pilot|medium|full)$", "", vapply(done, function(i) i$name, character(1)))
  chosen <- vapply(unique(family), function(f) {
    cand <- done[family == f]; nm <- vapply(cand, function(i) i$name, character(1))
    i <- .pick_mode(sub("^.*_", "", nm)); if (is.na(i)) NA_character_ else cand[[i]]$dir
  }, character(1))
  chosen <- chosen[!is.na(chosen)]
  if (!length(chosen)) stop("No finished ", Sys.getenv("SIM_MODE"), " run under ", root, call. = FALSE)
  skipped <- Filter(function(i) !i$finished, infos)
  if (length(skipped)) message("Not plotted:\n", paste(vapply(skipped, explain_folder, character(1)), collapse = "\n"))
  unname(chosen)
}

# ---- plotting one family ---------------------------------------------------------------------
plot_family <- function(results_dir) {
  if (!dir.exists(results_dir)) stop("Folder not found: ", results_dir, call. = FALSE)
  train_file <- list.files(results_dir, pattern = "^explore_.*_train\\.csv$", full.names = TRUE)
  if (length(train_file) != 1)
    stop("No finished explorer results in ", results_dir, ".\n", explain_folder(describe_folder(results_dir)), call. = FALSE)
  test_file <- sub("_train\\.csv$", "_test.csv", train_file)
  family <- sub("^explore_(.*)_train\\.csv$", "\\1", basename(train_file))
  mode_label <- sub("^.*_", "", basename(normalizePath(results_dir)))
  all_sets <- rbind(read.csv(train_file, stringsAsFactors = FALSE), read.csv(test_file, stringsAsFactors = FALSE))
  all_sets <- all_sets[all_sets$set %in% USE_SETS & all_sets$successful >= 20, ]
  if (nzchar(EXCLUDE_LABELS)) all_sets <- all_sets[!grepl(EXCLUDE_LABELS, all_sets$label), ]
  if ("p_rate1" %in% names(all_sets)) { all_sets$p_rate_average <- (all_sets$p_rate1 + all_sets$p_rate2) / 2; rate_var <- "p_rate_average" }
  else rate_var <- "p_rate"
  keys <- KEY_VARIABLES[[family]]
  keys <- keys[keys %in% names(all_sets)]
  if (length(keys) < 2) {
    numeric_cols <- setdiff(names(all_sets)[grepl("^[pd]_", names(all_sets)) & sapply(all_sets, is.numeric)], c(rate_var, "p_effect"))
    keys <- head(c(keys, numeric_cols), 2)
  }
  plot_dir <- file.path(results_dir, PLOT_FOLDER)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)        # created if missing
  if (!dir.exists(plot_dir)) stop("Cannot create ", plot_dir, call. = FALSE)

  open_device <- function(file, w = FIG_WIDTH_IN, h = FIG_HEIGHT_IN) {
    a <- list(filename = file, width = w, height = h, units = "in", res = FIG_DPI)
    if (FIG_FORMAT == "tiff") a$compression <- "lzw"
    if (FIG_FORMAT == "jpeg") a$quality <- 95
    if (isTRUE(capabilities("cairo"))) a$type <- "cairo"
    do.call(get(FIG_FORMAT), a)
  }
  safe <- function(x) gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", x))
  pretty <- function(v) gsub("_", " ", sub("^[pd]_", "", v))
  palette_for <- function(levels) setNames(c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")[seq_along(levels)], levels)
  # Splits a true-effect variable into three bands (or one band if it has too few distinct values).
  effect_groups <- function(x) {
    b <- unique(quantile(x, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE))
    if (length(b) < 3) return(factor(rep("all", length(x))))
    cut(x, b, include.lowest = TRUE, dig.lab = 2)
  }
  wide_range <- function(x) all(x > 0, na.rm = TRUE) && max(x, na.rm = TRUE) / min(x, na.rm = TRUE) > 20

  # Scatter of one outcome against one variable, coloured by mechanism, with a smooth per mechanism.
  scatter <- function(d, xv, yv, ylab, file, main, group_var = "p_mechanism", group_title = "mechanism") {
    if (nrow(d) < 20 || !xv %in% names(d)) return(invisible(FALSE))
    logx <- !grepl("rate", xv) && wide_range(d[[xv]]); mech <- factor(d[[group_var]]); pal <- palette_for(levels(mech))   # rates stay linear
    open_device(file); par(mar = c(4.5, 4.5, 3.5, 1))
    plot(d[[xv]], d[[yv]], log = if (logx) "x" else "", pch = 16, cex = 0.7, col = adjustcolor(pal[as.character(mech)], 0.6),
         xlab = pretty(xv), ylab = ylab, ylim = c(0, 1), main = main, cex.main = 0.95)
    if (yv == "rejection_rate") abline(h = c(0.05, HARM_LEVEL), lty = c(3, 2), col = "gray45")
    for (m in levels(mech)) {
      s <- d[mech == m & is.finite(d[[xv]]) & is.finite(d[[yv]]), ]
      if (nrow(s) >= 10) { x <- if (logx) log(s[[xv]]) else s[[xv]]; f <- lowess(x, s[[yv]], f = 0.6)
        lines(if (logx) exp(f$x) else f$x, f$y, col = pal[m], lwd = 2) }
    }
    legend("topright", levels(mech), col = pal, pch = 16, bty = "n", cex = 0.75, title = group_title)
    mtext(sprintf("%s, %s; %d configurations", family, mode_label, nrow(d)), side = 3, line = 0.2, cex = 0.7)
    dev.off(); invisible(TRUE)
  }
  # Average outcome in bins of two variables, one small map per mechanism (shared bins and colours).
  tile_map <- function(d, xv, yv, file, main) {
    if (nrow(d) < 40 || !all(c(xv, yv) %in% names(d))) return(invisible(FALSE))
    mechs <- sort(unique(d$p_mechanism)); k <- length(mechs)
    nb <- max(3L, min(6L, as.integer(floor(sqrt(nrow(d) / k / 6)))))
    edges <- function(v) unique(quantile(v, probs = seq(0, 1, length.out = nb + 1), na.rm = TRUE))
    bx <- edges(d[[xv]]); by <- edges(d[[yv]])
    if (length(bx) < 3 || length(by) < 3) return(invisible(FALSE))
    xi <- cut(d[[xv]], bx, include.lowest = TRUE, labels = FALSE); yi <- cut(d[[yv]], by, include.lowest = TRUE, labels = FALSE)
    lab <- function(e) format(signif(e, 2), scientific = FALSE, trim = TRUE)
    pal <- colorRampPalette(c("white", "#FDD49E", "#FC8D59", "#D7301F", "#7F0000"))(60)
    open_device(file, 10, 8); par(mfrow = c(ceiling(k / 2), min(k, 2)), mar = c(6.2, 5, 3, 1), oma = c(0, 0, 3, 0))
    for (mech in mechs) {
      m <- matrix(NA_real_, length(bx) - 1, length(by) - 1)
      for (i in seq_len(nrow(m))) for (j in seq_len(ncol(m))) {
        sel <- which(xi == i & yi == j & d$p_mechanism == mech)
        if (length(sel) >= 3) m[i, j] <- mean(d$rejection_rate[sel])
      }
      image(seq_len(nrow(m)), seq_len(ncol(m)), m, zlim = c(0, 0.6), axes = FALSE, xlab = "", ylab = "", col = pal,
            main = mech, cex.main = 0.95)
      axis(1, at = seq_len(nrow(m)), labels = paste0(lab(head(bx, -1)), "+"), cex.axis = 0.7, las = 2)
      axis(2, at = seq_len(ncol(m)), labels = paste0(lab(head(by, -1)), "+"), cex.axis = 0.7, las = 1)
      title(xlab = pretty(xv), line = 5); title(ylab = pretty(yv), line = 3.6)
      for (i in seq_len(nrow(m))) for (j in seq_len(ncol(m)))
        if (!is.na(m[i, j])) text(i, j, sprintf("%.2f", m[i, j]), cex = 0.65, col = if (m[i, j] > 0.35) "white" else "black")
      box()
    }
    mtext(main, outer = TRUE, line = 1.2, font = 2)
    mtext(sprintf("%s, %s; cell = mean Type I error (white = about 0.05 or lower; blank = too few configurations)", family, mode_label),
          outer = TRUE, line = -0.2, cex = 0.7)
    dev.off(); invisible(TRUE)
  }

  labels <- sort(unique(all_sets$label)); made <- 0L
  for (lab in labels) {
    d <- all_sets[all_sets$label == lab, ]; null <- d[d$p_effect == 0, ]; alt <- d[d$p_effect != 0, ]
    out <- file.path(plot_dir, safe(lab)); dir.create(out, recursive = TRUE, showWarnings = FALSE)
    made <- made + scatter(null, rate_var, "rejection_rate", "Type I error", file.path(out, sprintf("1_typeI_vs_missing_rate.%s", FIG_FORMAT)),
                           sprintf("%s: Type I error against the missing rate", lab))
    for (k in seq_along(keys))
      made <- made + scatter(null, keys[k], "rejection_rate", "Type I error", file.path(out, sprintf("%d_typeI_vs_%s.%s", k + 1L, safe(pretty(keys[k])), FIG_FORMAT)),
                             sprintf("%s: Type I error against %s", lab, pretty(keys[k])))
    made <- made + tile_map(null, keys[1], keys[2], file.path(out, sprintf("4_typeI_map_%s_by_%s.%s", safe(pretty(keys[1])), safe(pretty(keys[2])), FIG_FORMAT)),
                            sprintf("%s: where is the test safe?", lab))
    if (nrow(alt) >= 20) alt$effect_size <- effect_groups(alt$p_effect)
    made <- made + scatter(alt, keys[1], "rejection_rate", "Power", file.path(out, sprintf("5_power_vs_%s.%s", safe(pretty(keys[1])), FIG_FORMAT)),
                           sprintf("%s: power against %s, by size of the true effect", lab, pretty(keys[1])),
                           group_var = "effect_size", group_title = "true effect size")
    if (any(d$failure_rate > 0.01))
      made <- made + scatter(d, keys[1], "failure_rate", "Share of failed analyses", file.path(out, sprintf("6_failure_vs_%s.%s", safe(pretty(keys[1])), FIG_FORMAT)),
                             sprintf("%s: analyses that could not run, against %s", lab, pretty(keys[1])))
  }

  # ---- comparisons across the tests of the family ---------------------------------------------
  over <- file.path(plot_dir, "_overview"); dir.create(over, recursive = TRUE, showWarnings = FALSE)
  null_all <- all_sets[all_sets$p_effect == 0, ]
  if (nrow(null_all) >= 20) {
    open_device(file.path(over, sprintf("typeI_by_test_and_mechanism.%s", FIG_FORMAT)))
    grp <- interaction(null_all$label, null_all$p_mechanism, sep = " | ", drop = TRUE)
    longest <- max(strwidth(levels(grp), units = "inches", cex = 0.6))
    par(mar = c(longest / par("csi") + 2, 4.5, 3.5, 1))          # room for the longest label
    boxplot(null_all$rejection_rate ~ grp, las = 2, cex.axis = 0.6, ylim = c(0, 1), ylab = "Type I error", xlab = "",
            col = palette_for(sort(unique(null_all$p_mechanism)))[sub("^.* \\| ", "", levels(grp))], main = "Type I error by test and mechanism")
    abline(h = c(0.05, HARM_LEVEL), lty = c(3, 2), col = "gray45")
    mtext(sprintf("%s, %s; each box summarises random configurations", family, mode_label), side = 3, line = 0.2, cex = 0.7)
    dev.off(); made <- made + 1L
    open_device(file.path(over, sprintf("typeI_vs_missing_rate_all_tests.%s", FIG_FORMAT)))
    k <- length(labels); par(mfrow = c(ceiling(k / 3), min(k, 3)), mar = c(4, 4, 3, 1))
    pal <- palette_for(sort(unique(null_all$p_mechanism)))
    for (lab in labels) {
      d <- null_all[null_all$label == lab, ]
      plot(d[[rate_var]], d$rejection_rate, pch = 16, cex = 0.5, ylim = c(0, 1), col = adjustcolor(pal[as.character(d$p_mechanism)], 0.6),
           xlab = pretty(rate_var), ylab = "Type I error", main = lab, cex.main = 0.85)
      abline(h = c(0.05, HARM_LEVEL), lty = c(3, 2), col = "gray45")
    }
    dev.off(); made <- made + 1L
  }
  cat(sprintf("%s (%s): %d images in %s\n", family, mode_label, made, normalizePath(plot_dir)))
  invisible(made)
}

run_plots <- function() {
  if (!FIG_FORMAT %in% c("png", "tiff", "jpeg")) stop("FIG_FORMAT must be png, tiff or jpeg", call. = FALSE)
  targets <- if (nzchar(RESULTS_DIR)) RESULTS_DIR else find_runs(OUTPUT_ROOT)
  cat("Plotting:", paste(basename(targets), collapse = ", "), "\n\n")
  invisible(lapply(targets, plot_family))
}
# The whole workflow runs inside run_plots(), so a problem stops the run cleanly even if you paste
# line by line. Every folder is created if it is missing.
plot_results <- run_plots()
