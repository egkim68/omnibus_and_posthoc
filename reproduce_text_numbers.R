# REPRODUCE THE NUMBERS IN THE TEXT that the pipeline does not print. Base R only. Changes no result.
# Reads the per-design results of one ANOVA run and one post hoc run and prints, section by section:
#   Section 2.2  the missing-rate cap under group-dependent missingness and its effect on the rate-only AUC
#   Section 3.1  Welch ANOVA and Alexander-Govern by mechanism, smallest expected group below and above 4
#   Section 3.4  the frozen rule by slice (outliers, shape, missing rate, mechanism, selection strength)
#   Section 3.5  Table 2, the gated pipelines, the common subset, and the switching rule
#   Section 3.6  harmful shares at 0.075, 0.10 and 0.15, the post hoc slices, and the Dunn-Holm conditions
# Everything is also written to <output folder>/text_numbers_<mode>.txt.
# Use:  Rscript reproduce_text_numbers.R            (finds the runs of SIM_MODE, default medium)
#       Rscript reproduce_text_numbers.R <anova run folder> <post hoc run folder>

args <- commandArgs(trailingOnly = TRUE)
OUTPUT_ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "N:/paper-simulation-explorer")
MODE <- Sys.getenv("SIM_MODE", unset = "medium")
ANOVA_DIR <- if (length(args) >= 1) args[1] else file.path(OUTPUT_ROOT, paste0("anova_", MODE))
POSTHOC_DIR <- if (length(args) >= 2) args[2] else file.path(OUTPUT_ROOT, paste0("posthoc_", MODE))
SCRIPT_DIR <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); if (length(f)) dirname(normalizePath(f[1])) else file.path(OUTPUT_ROOT, "R") })
HARM <- 0.10
FROZEN <- "Selector BF<5 else Welch"

read_run <- function(d, fam) {
  a <- read.csv(file.path(d, sprintf("explore_%s_train.csv", fam)), stringsAsFactors = FALSE)
  b <- read.csv(file.path(d, sprintf("explore_%s_test.csv", fam)), stringsAsFactors = FALSE)
  a$set <- "train"; b$set <- "test"; x <- rbind(a, b); x <- x[x$successful >= 20, ]; x$key <- paste(x$set, x$config_id); x
}
pct <- function(x) sprintf("%.1f%%", 100 * x)
hs <- function(r, h = HARM) mean(r > h)
out_file <- file.path(OUTPUT_ROOT, sprintf("text_numbers_%s.txt", MODE))
sink(out_file, split = TRUE)
cat("Numbers in the text, from", normalizePath(ANOVA_DIR), "and", normalizePath(POSTHOC_DIR), "\n")

an <- read_run(ANOVA_DIR, "anova"); ph <- read_run(POSTHOC_DIR, "posthoc")
nul <- an[an$p_effect == 0, ]

# ---- Section 2.2: the 0.9 cap under group-dependent missingness ------------------------------------------
# The design generator of explore_anova.R is reused (only its make_design function is evaluated), so the
# group rates are rebuilt exactly from each design's seed.
cat("\n== Section 2.2: cap on group rates (group-dependent missingness, ANOVA designs) ==\n")
ex <- parse(file.path(SCRIPT_DIR, "explore_anova.R"))
is_md <- vapply(ex, function(e) is.call(e) && identical(as.character(e[[1]]), "<-") && identical(as.character(e[[2]]), "make_design"), logical(1))
eval(ex[[which(is_md)]])
des <- an[!duplicated(an$key), ]
mar <- des[des$p_mechanism == "MAR_group", ]
mar$eff_rate <- vapply(seq_len(nrow(mar)), function(i) { d <- make_design(mar[i, ]); sum(d$n_g * d$rate_g) / sum(d$n_g) }, numeric(1))
drop <- mar$p_rate - mar$eff_rate; bound <- drop > 1e-9
cat(sprintf("cap binds in %d of %d designs (%s); expected missing proportion lower by %.3f on average in those designs (%.3f over all of them)\n",
            sum(bound), nrow(mar), pct(mean(bound)), mean(drop[bound]), mean(drop)))
eff <- setNames(mar$eff_rate, mar$key)
auc <- function(score, y) { r <- rank(score); n1 <- sum(y); n0 <- sum(!y); if (!n1 || !n0) NA else (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
rate_auc <- function(x, rate) {
  tr <- x$set == "train"; y <- x$rejection_rate > HARM
  if (sum(y[tr]) < 5 || sum(y[!tr]) < 5) return(NA)
  fit <- suppressWarnings(glm(y ~ r + I(r^2), family = binomial, data = data.frame(y = y[tr], r = rate[tr])))
  auc(predict(fit, newdata = data.frame(r = rate[!tr])), y[!tr])
}
for (l in c("F test", "Welch ANOVA", "Brown-Forsythe", "Alexander-Govern", "Kruskal-Wallis")) {
  x <- nul[nul$label == l, ]; r2 <- ifelse(x$key %in% names(eff), eff[x$key], x$p_rate)
  cat(sprintf("  %-17s held-out rate-only AUC: target rate %.3f; expected rate after cap %.3f\n", l, rate_auc(x, x$p_rate), rate_auc(x, r2)))
}

# ---- Section 3.1 ------------------------------------------------------------------------------------------
cat("\n== Section 3.1: harmful share by mechanism, smallest expected group below 4 / 4 or more (null designs) ==\n")
for (l in c("Welch ANOVA", "Alexander-Govern")) for (m in c("MCAR", "MAR_group", "MNAR")) {
  x <- nul[nul$label == l & nul$p_mechanism == m, ]; s <- x$d_min_expected_n < 4
  cat(sprintf("  %-17s %-9s below 4: %s of %d (mean error %.3f); 4 or more: %s of %d\n", l, m, pct(hs(x$rejection_rate[s])), sum(s), mean(x$rejection_rate[s]), pct(hs(x$rejection_rate[!s])), sum(!s)))
}

# ---- Section 3.4 ------------------------------------------------------------------------------------------
cat("\n== Section 3.4: harmful share of the frozen rule, Welch ANOVA, F test and Brown-Forsythe by slice (null designs, pooled) ==\n")
slice_tab <- function(x, labels, by) {
  for (lev in unique(by[!is.na(by)])) {
    cells <- vapply(labels, function(l) { s <- x$label == l & by == lev & !is.na(by); sprintf("%s", pct(hs(x$rejection_rate[s]))) }, character(1))
    cat(sprintf("  %-22s n=%-4d %s\n", lev, sum(x$label == labels[1] & by == lev, na.rm = TRUE), paste(sprintf("%s %s", sub("^Selector ", "", labels), cells), collapse = "; ")))
  }
}
L4 <- c(FROZEN, "Welch ANOVA", "F test", "Brown-Forsythe")
slice_tab(nul, L4, ifelse(nul$p_outlier > 0, "outliers", "clean"))
slice_tab(nul, L4, nul$p_shape)
slice_tab(nul, L4, as.character(cut(nul$p_rate, c(-Inf, 0.1, 0.2, Inf), labels = c("rate <= 0.10", "rate 0.10 to 0.20", "rate > 0.20"))))
slice_tab(nul, L4, ifelse(nul$p_mechanism != "MNAR", nul$p_mechanism, ifelse(nul$p_strength <= 1.5, "MNAR weak (<=1.5)", "MNAR strong (>1.5)")))

# ---- Section 3.5 ------------------------------------------------------------------------------------------
cat("\n== Section 3.5 and Table 2 ==\n")
PROCS <- c("Tukey-Kramer", "Games-Howell", "Welch-Holm", "Dunn-Holm")
PIPES <- c("F > Tukey-Kramer", "Welch > Games-Howell", "Brown-Forsythe > Games-Howell", "Kruskal-Wallis > Dunn-Holm")
val <- function(prefix, p, col = "rejection_rate") { z <- ph[ph$label == paste(prefix, p), ]; setNames(z[[col]], z$key) }
pp_all <- lapply(setNames(c(PROCS, PIPES), c(PROCS, PIPES)), function(p) val("Pair power", p, "mean_estimate"))
fin <- Reduce(intersect, lapply(pp_all[PROCS], function(v) names(v)[is.finite(v)]))
for (p in c(PROCS, PIPES)) {
  f <- val("FWER", p); fp <- val("Partial FWER", p); ap <- val("Power any-pair", p); pp <- pp_all[[p]]
  cat(sprintf("  %-30s FWER %.3f, above 0.10 %s (n=%d) | false-pair %.3f | any-pair %.3f | per-pair %.3f (n=%d) | within 0.10: %d of %d\n",
              p, mean(f), pct(hs(f)), length(f), mean(fp), mean(ap), mean(pp[fin]), length(fin), sum(fp[fin] <= HARM), length(fin)))
}
cat(sprintf("  designs with a true difference: %d; with finite per-pair power for all procedures: %d\n", length(val("Power any-pair", "Welch-Holm")), length(fin)))
fpm <- sapply(PROCS, function(p) val("Partial FWER", p)[fin]); common <- fin[apply(fpm <= HARM, 1, all)]
cat(sprintf("  common subset (all four within 0.10): %d designs; per-pair power %s\n", length(common),
            paste(sprintf("%s %.3f", PROCS, sapply(PROCS, function(p) mean(pp_all[[p]][common]))), collapse = ", ")))
# switching rule: Games-Howell when the smallest EXPECTED group is 5 or more, Welch-Holm otherwise (design level)
mexp <- setNames(ph$d_min_expected_n, ph$key)
sw <- function(prefix, col = "rejection_rate", keys = NULL) { a <- val(prefix, "Welch-Holm", col); b <- val(prefix, "Games-Howell", col)
  k <- if (is.null(keys)) names(a) else keys; ifelse(mexp[k] >= 5, b[k], a[k]) }
cat(sprintf("  switching rule (Games-Howell when smallest expected group >= 5): FWER %.3f vs Welch-Holm %.3f; per-pair power %.3f vs %.3f\n",
            mean(sw("FWER")), mean(val("FWER", "Welch-Holm")), mean(sw("Pair power", "mean_estimate", fin)), mean(pp_all[["Welch-Holm"]][fin])))

# ---- Section 3.6 ------------------------------------------------------------------------------------------
cat("\n== Section 3.6: harmful shares at three lines ==\n")
for (h in c(0.075, 0.10, 0.15)) {
  cat(sprintf("  line %.3f  post hoc: %s\n", h, paste(sprintf("%s %s", PROCS, sapply(PROCS, function(p) pct(hs(val("FWER", p), h)))), collapse = "; ")))
  cat(sprintf("             omnibus:  %s\n", paste(sprintf("%s %s", sub("^Selector ", "", L4[1:2]), sapply(L4[1:2], function(l) pct(hs(nul$rejection_rate[nul$label == l], h)))), collapse = "; ")))
}
cat("\n== Section 3.6: post hoc harmful share by slice (null designs, pooled) ==\n")
pn <- ph[ph$p_effect == 0 & ph$label %in% paste("FWER", PROCS), ]; pn$label <- sub("^FWER ", "", pn$label)
slice_tab(pn, PROCS, pn$p_mechanism)
slice_tab(pn, PROCS, pn$p_shape)
slice_tab(pn, PROCS, ifelse(pn$p_sd_ratio > 1, "unequal spread", "equal spread"))
slice_tab(pn, PROCS, ifelse(pn$p_mechanism != "MNAR", pn$p_mechanism, ifelse(pn$p_strength <= 1.5, "MNAR weak (<=1.5)", "MNAR strong (>1.5)")))
slice_tab(pn, PROCS, ifelse(pn$p_outlier == 0, "no outliers", ifelse(pn$p_outlier_size < 8, "outliers 4 to <8 SD", "outliers 8 to 12 SD")))
slice_tab(pn, PROCS, as.character(cut(pn$p_rate, c(-Inf, 0.1, 0.2, Inf), labels = c("rate <= 0.10", "rate 0.10 to 0.20", "rate > 0.20"))))

cat("\n== Section 3.5: Dunn-Holm and Welch-Holm, equal means, skewed data, unequal spread (mean familywise error) ==\n")
for (o in c(FALSE, TRUE)) {
  x <- pn[pn$p_shape == "skewed" & pn$p_sd_ratio > 1 & (pn$p_outlier > 0) == o, ]
  cat(sprintf("  %-12s Dunn-Holm %.3f, Welch-Holm %.3f (%d designs)\n", if (o) "outliers" else "no outliers",
              mean(x$rejection_rate[x$label == "Dunn-Holm"]), mean(x$rejection_rate[x$label == "Welch-Holm"]), sum(x$label == "Dunn-Holm")))
}
cat("\n== Section 3.5: per-pair power of Dunn-Holm and Welch-Holm by data shape, with and without outliers ==\n")
pa <- ph[ph$p_effect > 0 & ph$label %in% c("Pair power Dunn-Holm", "Pair power Welch-Holm"), ]
for (o in c(FALSE, TRUE)) for (sh in c("normal", "heavy", "skewed")) {
  x <- pa[(pa$p_outlier > 0) == o & pa$p_shape == sh, ]
  cat(sprintf("  %-12s %-7s Dunn-Holm %.2f, Welch-Holm %.2f\n", if (o) "outliers" else "no outliers", sh,
              mean(x$mean_estimate[x$label == "Pair power Dunn-Holm"], na.rm = TRUE), mean(x$mean_estimate[x$label == "Pair power Welch-Holm"], na.rm = TRUE)))
}
cat("\nAUCs at other harmful lines: rerun analyze_explorer.R with SIM_HARM_LEVEL set to 0.075 or 0.15 (into a copy of the run folder).\n")
sink()
cat("\nWritten to", normalizePath(out_file), "\n")
