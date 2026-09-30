# FIGURE 5 for the paper: the post hoc procedures. Base R only. Reads one post hoc explorer run
# (explore_posthoc_train.csv, explore_posthoc_test.csv) and writes ONE image, Figure_5_posthoc.png, to
# <run>/paper_figures. Companion to make_paper_figures.R (same settings, same visual style); run this
# after that script, or on its own. Rerun any time, e.g. after the full run.
#
# Figure 5 has two panels:
#   A) familywise error (no true difference) of the four post hoc procedures against the smallest
#      observed group, one line per procedure, pooling TRAIN and TEST for a stable picture.
#   B) validity-conditional per-pair power: average power among the designs (with a true difference in
#      some groups) where a procedure's own false-pair rate is at or below HARM_LEVEL, next to its raw
#      (not size-adjusted) power over all eligible designs, for TEST designs only (the held-out check).
#      The kept subset differs between procedures, so the gap between the two bars reflects both the
#      procedure and which designs were kept; it is not a comparison on identical designs.

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
FIG_FOLDER <- "paper_figures"
FIG_FORMAT <- "png"
FIG_WIDTH_IN <- 7.2
FIG_DPI <- 300
HARM_LEVEL <- 0.10
PROCS <- c("Tukey-Kramer", "Games-Howell", "Welch-Holm", "Dunn-Holm")
PROC_COLORS <- c(`Tukey-Kramer` = "#D55E00", `Games-Howell` = "#E69F00", `Welch-Holm` = "#0072B2", `Dunn-Holm` = "#009E73")
PROC_PCH <- c(`Tukey-Kramer` = 24, `Games-Howell` = 22, `Welch-Holm` = 21, `Dunn-Holm` = 23)

find_posthoc_run <- function(root) {
  if (!dir.exists(root)) stop("Output folder not found: ", root, call. = FALSE)
  dirs <- list.dirs(root, recursive = FALSE)
  dirs <- dirs[grepl("^posthoc_(pilot|medium|full)$", basename(dirs))]
  ok <- Filter(function(d) { f <- file.path(d, "explore_posthoc_train.csv"); g <- file.path(d, "explore_posthoc_test.csv"); file.exists(f) && file.exists(g) }, dirs)
  if (!length(ok)) stop("No finished post hoc run under ", root, ". Run explore_posthoc.R first.", call. = FALSE)
  nm <- sub("^posthoc_", "", basename(ok))
  i <- .pick_mode(nm)
  if (is.na(i)) stop("No finished post hoc ", Sys.getenv("SIM_MODE"), " run under ", root, call. = FALSE)
  ok[[i]]
}
panel_letter <- function(letter, line = 0.3) mtext(letter, side = 3, adj = 0, line = line, font = 2, cex = 0.95)

make_figure5 <- function(results_dir) {
  tr <- read.csv(file.path(results_dir, "explore_posthoc_train.csv"), stringsAsFactors = FALSE)
  te <- read.csv(file.path(results_dir, "explore_posthoc_test.csv"), stringsAsFactors = FALSE)
  tr$set <- "train"; te$set <- "test"
  d <- rbind(tr, te); d <- d[d$successful >= 20, ]; d$key <- paste(d$set, d$config_id)
  out_dir <- file.path(results_dir, FIG_FOLDER); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  f <- file.path(out_dir, sprintf("Figure_5_posthoc.%s", FIG_FORMAT))
  a <- list(filename = f, width = FIG_WIDTH_IN, height = 4.3, units = "in", res = FIG_DPI)
  if (isTRUE(capabilities("cairo"))) a$type <- "cairo"
  do.call(png, a)
  layout(matrix(1:2, 1, 2), widths = c(1.05, 1))

  # ---- Panel A: familywise error against smallest observed group, null designs, train + test --------
  par(mar = c(5.0, 3.8, 1.4, 0.6), mgp = c(2.3, 0.7, 0), tcl = -0.3, cex.axis = 1.05, cex.lab = 1.1, las = 1)
  null_d <- d[d$p_effect == 0 & grepl("^FWER ", d$label) & sub("^FWER ", "", d$label) %in% PROCS, ]
  null_d$proc <- sub("^FWER ", "", null_d$label)
  plot(NA, xlim = c(0.3, 90), ylim = c(0, 0.5), log = "x", xlab = "", ylab = "Familywise error rate", axes = FALSE)
  rect(0.1, -0.05, 200, HARM_LEVEL, col = rgb(0, 0, 0, 0.05), border = NA)
  abline(h = 0.05, lty = 3, col = "gray40"); abline(h = HARM_LEVEL, lty = 2, col = "gray40")
  axis(1, at = c(0.5, 1, 2, 5, 10, 20, 50), labels = c("0.5", "1", "2", "5", "10", "20", "50"), cex.axis = 0.95)
  title(xlab = "Smallest expected group\nafter missingness (log scale)", line = 3.4); axis(2, at = seq(0, 0.5, 0.1), cex.axis = 1.05); box()
  for (p in PROCS) {
    x <- null_d[null_d$proc == p, ]
    if (nrow(x) >= 10) { fit <- lowess(log(pmax(x$d_min_expected_n, 0.3)), x$rejection_rate, f = 0.6); lines(exp(fit$x), fit$y, col = PROC_COLORS[[p]], lwd = 2.2) }
  }
  panel_letter("A", 0.4)
  legend("topright", legend = c(PROCS, "Nominal 0.05", sprintf("Harmful above %.2f", HARM_LEVEL)),
         col = c(PROC_COLORS[PROCS], "gray40", "gray40"), lty = c(rep(1, 4), 3, 2), lwd = c(rep(2.4, 4), 1.2, 1.2), bty = "n", cex = 0.95, seg.len = 1.8)

  # ---- Panel B: validity-conditional vs raw per-pair power, held-out (test) designs with a true difference --
  par(mar = c(5.0, 3.8, 1.4, 0.6), cex.axis = 1.05, cex.lab = 1.1)
  alt <- d[d$set == "test" & d$p_effect > 0, ]
  get <- function(prefix, col) { x <- alt[alt$label %in% paste(prefix, PROCS), c("key", "label", col)]; x$proc <- sub(paste0("^", prefix, " "), "", x$label)
    setNames(lapply(PROCS, function(p) setNames(x[[col]][x$proc == p], x$key[x$proc == p])), PROCS) }
  pf <- get("Partial FWER", "rejection_rate"); pp <- get("Pair power", "mean_estimate")
  raw_power <- sapply(PROCS, function(p) mean(pp[[p]], na.rm = TRUE))
  cond_power <- sapply(PROCS, function(p) { keys <- names(pf[[p]])[!is.na(pf[[p]]) & pf[[p]] <= HARM_LEVEL]; mean(pp[[p]][keys], na.rm = TRUE) })
  n_kept <- sapply(PROCS, function(p) sum(!is.na(pf[[p]]) & pf[[p]] <= HARM_LEVEL))
  mat <- rbind(raw = raw_power, valid = cond_power)
  ymax <- max(mat, na.rm = TRUE) * 1.75
  plot(NA, xlim = c(0.3, length(PROCS) + 0.7), ylim = c(0, ymax), axes = FALSE, xlab = "", ylab = "Per-pair power")
  axis(2, cex.axis = 1.05); box()
  bw <- 0.32
  for (i in seq_along(PROCS)) {
    p <- PROCS[i]; x0 <- i
    rect(x0 - bw, 0, x0, raw_power[p], col = "gray75", border = "white")
    rect(x0, 0, x0 + bw, cond_power[p], col = PROC_COLORS[[p]], border = "white")
    text(x0 + bw / 2, cond_power[p] + 0.05 * ymax, sprintf("n=%d", n_kept[p]), cex = 0.85, col = "gray30")
    text(x0, -0.04 * ymax, sub("-", "-\n", p), xpd = TRUE, cex = 0.85, adj = c(0.5, 1))
  }
  panel_letter("B", 0.4)
  # the second key is a strip of the four procedure colors, matching the colored bars
  lg <- legend("topleft", legend = c("Raw power (not size-adjusted)", sprintf("Validity-conditional power\n(error \u2264 %.2f)", HARM_LEVEL)),
               fill = c("gray75", NA), border = NA, bty = "n", cex = 0.95, y.intersp = 1.6)
  kw <- (lg$text$x[2] - lg$rect$left) * 0.45; kh <- 0.022 * ymax; kx <- lg$text$x[2] - kw * 1.35
  for (j in seq_along(PROCS)) rect(kx + (j - 1) * kw / 4, lg$text$y[2] - kh, kx + j * kw / 4, lg$text$y[2] + kh, col = PROC_COLORS[[PROCS[j]]], border = NA)

  dev.off()
  cat(sprintf("Figure 5 written to %s\n", normalizePath(f)))
  invisible(f)
}

run_figure5 <- function() {
  target <- if (nzchar(RESULTS_DIR)) RESULTS_DIR else find_posthoc_run(OUTPUT_ROOT)
  cat("Building Figure 5 from:", basename(target), "\n")
  make_figure5(target)
}
figure5_path <- run_figure5()
