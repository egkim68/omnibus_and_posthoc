# CHECK THE PRE-STATED HYPOTHESES (Appendix B of the manuscript; PROTOCOL_v1)
# Base R only. New in version 1.5. Reads the finished runs of one mode (ANOVA and/or post hoc) and their
# "analysis" folders, and writes one table per family with every criterion, the numbers behind it and a
# verdict:  <run>/analysis/confirmation_omnibus.csv  and  <run>/analysis/confirmation_posthoc.csv.
# It changes no result and runs no simulation: rerun it any time. In a medium run the table is only a
# preview (the exploration generated the hypotheses, so it cannot test them); in the full run it is the test.
#
# HOW EACH CRITERION IS READ (so the reading is on record, not decided after seeing the numbers):
#   H1   per procedure, on the held-out designs (analysis_validation.csv): the AUC of the missing rate alone
#        is below 0.70 AND the design-aware model's AUC is at least 0.15 higher. The protocol requires this
#        for at least 5 of its 6 key procedures; the table reports every procedure and the count.
#   H2   null designs, train and test pooled, cells as in Figure 2 (precision channel: smallest expected
#        group below 4; bias channel: MNAR): Welch's precision-only harmful share exceeds Brown-Forsythe's
#        by at least 20 points; Brown-Forsythe's bias-only share is at least 15% and at least 10 points
#        above its precision-only share.
#   H3   null designs, train and test pooled ("all full designs"): the FROZEN rule's harmful share is at
#        most 15% and at least 10 points below both Welch ANOVA's and the F test's.
#   H4   null designs, train and test pooled: Welch-Holm's familywise-error harmful share is at most 15%
#        and at least 10 points below both Games-Howell's and Tukey-Kramer's.
#   H4b  (a) designs with outliers and a true difference, pooled: Dunn-Holm's mean per-pair power exceeds
#        Welch-Holm's by at least 0.05, and Dunn-Holm's harmful share in null designs with outliers is at
#        most 15%; (b) null designs with skewed data and unequal spread (SD ratio above 1): Dunn-Holm's
#        harmful share exceeds Welch-Holm's.
# "Harmful" means an error rate above HARM_LEVEL (0.10). Only configurations with at least 20 successful
# replications are used, as in the other scripts.

OUTPUT_ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "N:/paper-simulation-explorer")
CHECK_MODE <- Sys.getenv("SIM_MODE", unset = "full")
HARM_LEVEL <- 0.10
PRECISION_LIMIT <- 4
# H1 key procedures as named in docs/PROTOCOL_v1.md (labels as in analysis_validation.csv). If left empty,
# the verdict is given for the worst case: it holds for ANY six only if at most one procedure fails or
# cannot be estimated.
H1_KEY_PROCEDURES <- character(0)
# every procedure row H1 is scored on (5 omnibus tests, 4 post hoc procedures, 4 gated pipelines)
H1_EXPECTED <- c("Alexander-Govern", "Brown-Forsythe", "F test", "Kruskal-Wallis", "Welch ANOVA",
                 paste("FWER", c("Tukey-Kramer", "Games-Howell", "Welch-Holm", "Dunn-Holm", "F > Tukey-Kramer",
                                 "Welch > Games-Howell", "Brown-Forsythe > Games-Howell", "Kruskal-Wallis > Dunn-Holm")))
FROZEN_RULE <- "Selector BF<5 else Welch"   # frozen before the confirmatory run (PROTOCOL_v1); do not change

check_hypotheses <- function(root, mode) {
  if (!dir.exists(root)) stop("Output folder not found: ", root, call. = FALSE)
  read_run <- function(family) {
    d <- file.path(root, paste0(family, "_", mode))
    tr <- file.path(d, sprintf("explore_%s_train.csv", family)); te <- file.path(d, sprintf("explore_%s_test.csv", family))
    if (!file.exists(tr) || !file.exists(te)) return(NULL)
    a <- read.csv(tr, stringsAsFactors = FALSE); b <- read.csv(te, stringsAsFactors = FALSE); a$set <- "train"; b$set <- "test"
    x <- rbind(a, b); x <- x[x$successful >= 20, ]
    list(dir = d, data = x, validation = { v <- file.path(d, "analysis", "analysis_validation.csv"); if (file.exists(v)) read.csv(v, stringsAsFactors = FALSE) else NULL })
  }
  harm <- function(x, label) { r <- x$rejection_rate[x$label == label]; c(k = sum(r > HARM_LEVEL, na.rm = TRUE), n = sum(!is.na(r))) }
  share <- function(h) if (h[["n"]] > 0) h[["k"]] / h[["n"]] else NA_real_
  row <- function(h, part, value, criterion, pass) data.frame(hypothesis = h, part = part, value = value, criterion = criterion,
    verdict = if (is.na(pass)) "not estimable" else if (pass) "MET" else "NOT MET", stringsAsFactors = FALSE)
  pct <- function(p) sprintf("%.1f%%", 100 * p)
  h1_rows <- function(v) {
    if (is.null(v)) return(NULL)
    v <- v[v$truth %in% "null", ]
    do.call(rbind, lapply(seq_len(nrow(v)), function(i) {
      r <- v$auc_rate_only[i]; m <- v$auc_model[i]
      row("H1", v$label[i], sprintf("missing-rate AUC %.3f; model AUC %.3f; gain %.3f", r, m, m - r), "rate AUC < 0.70 and gain >= 0.15",
          if (is.na(r) || is.na(m)) NA else (r < 0.70 && m - r >= 0.15))
    }))
  }
  finish <- function(tab, run, file) {
    stage <- if (mode == "full") "CONFIRMATORY TEST" else "PREVIEW ONLY (exploratory data; not a test)"
    tab$stage <- stage
    out <- file.path(run$dir, "analysis"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
    write.csv(tab, file.path(out, file), row.names = FALSE)
    cat(sprintf("\n%s  [%s]  written to %s\n", file, stage, normalizePath(file.path(out, file))))
    for (i in seq_len(nrow(tab))) cat(sprintf("  %-4s %-34s %-9s %s\n", tab$hypothesis[i], substr(tab$part[i], 1, 34), tab$verdict[i], tab$value[i]))
    tab
  }
  results <- list()

  an <- read_run("anova")
  if (!is.null(an)) {
    x <- an$data; nul <- x[x$p_effect == 0, ]
    tab <- h1_rows(an$validation)
    # H2
    cell <- function(l, prec, bias) { s <- nul[nul$label == l & (nul$d_min_expected_n < PRECISION_LIMIT) == prec & (nul$p_mechanism == "MNAR") == bias, ]
      c(k = sum(s$rejection_rate > HARM_LEVEL), n = nrow(s)) }
    w_p <- cell("Welch ANOVA", TRUE, FALSE); b_p <- cell("Brown-Forsythe", TRUE, FALSE); b_b <- cell("Brown-Forsythe", FALSE, TRUE)
    g1 <- share(w_p) - share(b_p)
    tab <- rbind(tab,
      row("H2", "Welch vs BF, precision only", sprintf("Welch %s (%d/%d); BF %s (%d/%d); gap %.1f points", pct(share(w_p)), w_p[["k"]], w_p[["n"]], pct(share(b_p)), b_p[["k"]], b_p[["n"]], 100 * g1),
          "gap >= 20 points", g1 >= 0.20),
      row("H2", "BF bias only vs precision only", sprintf("bias only %s (%d/%d); precision only %s; gap %.1f points", pct(share(b_b)), b_b[["k"]], b_b[["n"]], pct(share(b_p)), 100 * (share(b_b) - share(b_p))),
          "bias-only share >= 15% and >= 10 points above precision-only", share(b_b) >= 0.15 && share(b_b) - share(b_p) >= 0.10))
    # H3
    if (FROZEN_RULE %in% nul$label) {
      r <- harm(nul, FROZEN_RULE); w <- harm(nul, "Welch ANOVA"); f <- harm(nul, "F test")
      tab <- rbind(tab, row("H3", paste("frozen rule", sub("^Selector ", "", FROZEN_RULE)),
        sprintf("rule %s (%d/%d); Welch %s; F test %s; gaps %.1f and %.1f points", pct(share(r)), r[["k"]], r[["n"]], pct(share(w)), pct(share(f)), 100 * (share(w) - share(r)), 100 * (share(f) - share(r))),
        "rule share <= 15% and >= 10 points below Welch and F test", share(r) <= 0.15 && share(w) - share(r) >= 0.10 && share(f) - share(r) >= 0.10))
      te <- nul[nul$set == "test", ]; rt <- harm(te, FROZEN_RULE); wt <- harm(te, "Welch ANOVA"); ft <- harm(te, "F test")
      tab <- rbind(tab, row("H3 (held-out only, descriptive)", "same rule, test designs only",
        sprintf("rule %s (%d/%d); Welch %s; F test %s", pct(share(rt)), rt[["k"]], rt[["n"]], pct(share(wt)), pct(share(ft))),
        "not the pre-stated criterion; reported for context", NA))
      tab$verdict[nrow(tab)] <- "descriptive"
    }
    results$omnibus <- finish(tab, an, "confirmation_omnibus.csv")
  }

  ph <- read_run("posthoc")
  if (!is.null(ph)) {
    x <- ph$data; nul <- x[x$p_effect == 0, ]; alt <- x[x$p_effect != 0, ]
    tab <- h1_rows(ph$validation)
    wh <- harm(nul, "FWER Welch-Holm"); gh <- harm(nul, "FWER Games-Howell"); tk <- harm(nul, "FWER Tukey-Kramer")
    tab <- rbind(tab, row("H4", "Welch-Holm vs Games-Howell, Tukey-Kramer",
      sprintf("Welch-Holm %s (%d/%d); Games-Howell %s; Tukey-Kramer %s", pct(share(wh)), wh[["k"]], wh[["n"]], pct(share(gh)), pct(share(tk))),
      "Welch-Holm <= 15% and >= 10 points below both", share(wh) <= 0.15 && share(gh) - share(wh) >= 0.10 && share(tk) - share(wh) >= 0.10))
    ao <- alt[alt$p_outlier > 0, ]
    # Per-pair power is conditional on a truly different pair surviving deletion in every replication of a
    # design (otherwise the design's mean pair power is NA). Both procedures are compared on the SAME designs,
    # those where both power estimates are finite, and that number is reported next to the eligible number.
    key_of <- function(p) { z <- ao[ao$label == paste("Pair power", p), ]; setNames(z$mean_estimate, paste(z$set, z$config_id)) }
    pd <- key_of("Dunn-Holm"); pwh <- key_of("Welch-Holm")
    both <- intersect(names(pd)[is.finite(pd)], names(pwh)[is.finite(pwh)])
    n_elig <- length(unique(paste(ao$set, ao$config_id)))
    pw <- function(p) mean((if (p == "Dunn-Holm") pd else pwh)[both])
    dn_o <- harm(nul[nul$p_outlier > 0, ], "FWER Dunn-Holm")
    gain <- pw("Dunn-Holm") - pw("Welch-Holm")
    tab <- rbind(tab, row("H4b", "(a) with outliers: power gain, Dunn safety",
      sprintf("per-pair power Dunn %.3f vs Welch-Holm %.3f (gain %.3f, %d of %d eligible designs with finite power); Dunn harmful with outliers %s (%d/%d)",
              pw("Dunn-Holm"), pw("Welch-Holm"), gain, length(both), n_elig, pct(share(dn_o)), dn_o[["k"]], dn_o[["n"]]),
      "gain >= 0.05 and Dunn harmful share <= 15%", gain >= 0.05 && share(dn_o) <= 0.15))
    sk <- nul[nul$p_shape == "skewed" & nul$p_sd_ratio > 1, ]
    ds <- harm(sk, "FWER Dunn-Holm"); ws <- harm(sk, "FWER Welch-Holm")
    tab <- rbind(tab, row("H4b", "(b) equal means, skewed, unequal spread",
      sprintf("Dunn %s (%d/%d) vs Welch-Holm %s (%d/%d)", pct(share(ds)), ds[["k"]], ds[["n"]], pct(share(ws)), ws[["k"]], ws[["n"]]),
      "Dunn harmful share > Welch-Holm's", share(ds) > share(ws)))
    results$posthoc <- finish(tab, ph, "confirmation_posthoc.csv")
  }
  if (!length(results)) cat("No finished", mode, "runs to check under", root, "\n")
  h1 <- do.call(rbind, lapply(results, function(t) t[t$hypothesis == "H1", ]))
  if (!is.null(h1) && nrow(h1)) {
    # H1 spans both families: a verdict is given only when every expected procedure row is present
    missing_rows <- setdiff(H1_EXPECTED, h1$part)
    n_met <- sum(h1$verdict == "MET"); n_na <- sum(h1$verdict == "not estimable"); bad <- sum(h1$verdict != "MET")
    if (length(missing_rows)) {
      overall <- "INCOMPLETE"; basis <- paste("missing rows:", paste(missing_rows, collapse = ", "))
    } else if (length(H1_KEY_PROCEDURES)) {
      k <- h1[h1$part %in% H1_KEY_PROCEDURES, ]
      if (nrow(k) != length(H1_KEY_PROCEDURES)) { overall <- "INCOMPLETE"; basis <- paste("key procedures not found:", paste(setdiff(H1_KEY_PROCEDURES, k$part), collapse = ", ")) }
      else { overall <- if (sum(k$verdict == "MET") >= 5) "MET" else "NOT MET"; basis <- sprintf("%d of %d key procedures met", sum(k$verdict == "MET"), nrow(k)) }
    } else {
      overall <- if (bad <= 1) "MET" else "UNDETERMINED (depends on which six were named)"
      basis <- sprintf("any six of %d procedures: %d not met or not estimable, at most 1 allowed", nrow(h1), bad)
    }
    cat(sprintf("\nH1 overall: %s. %d of %d procedures met the per-procedure criterion (%d not estimable); %s.\n", overall, n_met, nrow(h1), n_na, basis))
    h1_out <- data.frame(hypothesis = "H1 overall", met = n_met, rows = nrow(h1), not_estimable = n_na, basis = basis, verdict = overall,
                         stage = if (mode == "full") "CONFIRMATORY TEST" else "PREVIEW ONLY (exploratory data; not a test)", stringsAsFactors = FALSE)
    f_h1 <- file.path(root, sprintf("confirmation_H1_overall_%s.csv", mode))
    write.csv(h1_out, f_h1, row.names = FALSE); cat("  written to", normalizePath(f_h1), "\n")
    results$h1_overall <- h1_out
  }
  invisible(results)
}
hypothesis_results <- check_hypotheses(OUTPUT_ROOT, CHECK_MODE)
