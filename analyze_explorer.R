# EXPLORER ANALYSIS: what drives failures, which rules predict them, do the rules hold up?
# Works with the output of ANY of the four explorers (correlation, twogroup, anova, chisq).
# Base R plus the rpart package that ships with standard R installations.
# Reads explore_<family>_train.csv and explore_<family>_test.csv from a results folder and
# writes tables (CSV) and images (PNG) to an "analysis" folder inside it (created if missing).
# Independent of the simulations: rerun any time, for example with a different HARM_LEVEL.
#
# Use in R or RStudio: just source() or paste this file. It looks under OUTPUT_ROOT for every
# FINISHED explorer run and analyses each one (a "full" run is preferred over a "pilot" run of
# the same family). To analyse a single folder instead, set RESULTS_DIR below (or the
# environment variable SIM_RESULTS_DIR), or from a terminal:
#   Rscript analyze_explorer.R "N:/.../anova_full"
# If nothing can be analysed, the message lists what it found in each folder and why.
#
# WHAT IT DOES, for each analysis label (each test or handling), on NULL configurations
# (true effect zero, so rejection is a Type I error):
#  1 DRIVERS   a quasi-binomial model of the rejection rate on the configuration parameters;
#              each term's importance is the share of the model's null deviance that would
#              be lost by dropping it. The four strongest main effects are then crossed with one
#              another to look for interactions.
#  2 RULES     a small decision tree (rpart, depth at most 4) finds regions of the design space
#              where the Type I error exceeds HARM_LEVEL. Each region is written out as a
#              plain rule, with how often it held on the TRAIN set and on the fresh TEST set.
#  3 CHECK     on the TEST set it compares three ways of predicting which configurations are
#              harmful: the missing rate alone, the fitted model, and the tree. The area under
#              the ROC curve (AUC; 0.5 = no better than chance, 1 = perfect) summarises each.
#              If the missing rate alone already predicts harm as well as the others, the
#              percentage is a good guide; if not, layout and mechanism matter more.
#  4 OBSERVABLE  the same tree fitted only on what an analyst can observe (group sizes after
#              missingness, missing rates, spread, shape, outliers; NOT the mechanism, its strength
#              or its direction). The gap between its AUC and the full tree's AUC is the part of the
#              danger that cannot be seen from the data.
# AUCs and rules are skipped for a test with too few harmful configurations (see MIN_HARMFUL_*),
# because with a handful of cases they are noise; the validation table says so in its note column.
# It also lists the drivers of POWER on the configurations with a true effect.
# Findings are exploratory until they replicate on the TEST set; the TEST set is never used to
# choose rules, only to check them.

args <- commandArgs(trailingOnly = TRUE)
# Leave RESULTS_DIR empty to analyse every finished run found under OUTPUT_ROOT.
RESULTS_DIR <- if (length(args) >= 1) args[1] else Sys.getenv("SIM_RESULTS_DIR", unset = "")

# Which run to use: the mode set by the pipeline (SIM_MODE) if any, otherwise full > medium > pilot.
.pick_mode <- function(modes) {
  m <- Sys.getenv("SIM_MODE", unset = "")
  pref <- if (m %in% c("pilot", "medium", "full")) m else c("full", "medium", "pilot")
  i <- which.min(match(modes, pref))
  if (length(i)) i else NA_integer_
}
OUTPUT_ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "N:/paper-simulation-explorer")
HARM_LEVEL <- as.numeric(Sys.getenv("SIM_HARM_LEVEL", unset = "0.10"))   # a null configuration is "harmful" if its Type I error exceeds this
MIN_CONFIGS <- 40L        # fewest null configurations needed to analyse a label
TOP_TERMS <- 4L           # strongest main effects crossed to look for interactions
# Parameters an analyst cannot see in real data (how the values went missing). A second tree is
# fitted WITHOUT them, to show how much of the danger can be recognised from observable information.
UNOBSERVABLE <- c("p_mechanism", "p_strength", "p_direction", "p_alt_share")
# Labels matching this pattern are left out of the per-test analysis: the ANOVA "Selector ..." rules
# (evaluate_selectors.R handles them) and the post hoc power-side measures ("Partial FWER", "Power
# any-pair", "Pair power"), which have no Type I error reading. The post hoc "FWER ..." labels stay.
EXCLUDE_LABELS <- "^(Selector|Partial FWER|Power any-pair|Pair power)"
MIN_HARMFUL_TRAIN <- 10L  # AUCs and tree rules need at least this many harmful TRAIN configurations
MIN_HARMFUL_TEST <- 5L    # ... and at least this many harmful TEST configurations
# Output goes to this folder inside each results folder. Each run OVERWRITES the files in it, so
# to keep a version (for example before changing HARM_LEVEL) use a different name next time,
# such as "analysis_harm05".
ANALYSIS_FOLDER <- "analysis"
FIG_FORMAT <- "png"       # "png", "tiff" or "jpeg"
FIG_WIDTH_IN <- 11
FIG_HEIGHT_IN <- 7
FIG_DPI <- 300

# What is in an explorer folder: finished runs have exactly one train and one test CSV.
describe_folder <- function(d) {
  train <- list.files(d, pattern = "^explore_.*_train\\.csv$")
  test <- list.files(d, pattern = "^explore_.*_test\\.csv$")
  checkpoints <- length(list.files(file.path(d, "checkpoints")))
  list(dir = d, name = basename(d), finished = length(train) == 1 && length(test) == 1,
       n_train = length(train), n_test = length(test), checkpoints = checkpoints)
}
explain_folder <- function(info) {
  what <- if (info$finished) "finished" else if (info$checkpoints > 0)
    sprintf("UNFINISHED: %d checkpoints but no results CSV yet; paste the explorer script again with the same settings to resume", info$checkpoints)
  else if (info$n_train + info$n_test > 0) "incomplete results files"
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
  if (!length(infos))
    stop("No explorer folders found under ", root, ". Run an explorer script first.", call. = FALSE)
  done <- Filter(function(i) i$finished, infos)
  if (!length(done))
    stop("No finished explorer run under ", root, ". What is there:\n",
         paste(vapply(infos, explain_folder, character(1)), collapse = "\n"), call. = FALSE)
  family <- sub("_(pilot|medium|full)$", "", vapply(done, function(i) i$name, character(1)))
  chosen <- vapply(unique(family), function(f) {
    cand <- done[family == f]; names_c <- vapply(cand, function(i) i$name, character(1))
    i <- .pick_mode(sub("^.*_", "", names_c)); if (is.na(i)) NA_character_ else cand[[i]]$dir
  }, character(1))
  chosen <- chosen[!is.na(chosen)]
  if (!length(chosen)) stop("No finished ", Sys.getenv("SIM_MODE"), " run under ", root, call. = FALSE)
  skipped <- Filter(function(i) !i$finished, infos)
  if (length(skipped)) message("Not analysed:\n", paste(vapply(skipped, explain_folder, character(1)), collapse = "\n"))
  unname(chosen)
}

analyze_explorer <- function(results_dir) {
  if (!FIG_FORMAT %in% c("png", "tiff", "jpeg"))
    stop("FIG_FORMAT must be png, tiff or jpeg", call. = FALSE)
  if (!requireNamespace("rpart", quietly = TRUE))
    stop("The rpart package is needed (it ships with standard R installations).", call. = FALSE)
  if (!dir.exists(results_dir))
    stop("Folder not found: ", results_dir, ". Check the path, or leave RESULTS_DIR empty to ",
         "find finished runs automatically.", call. = FALSE)
  train_file <- list.files(results_dir, pattern = "^explore_.*_train\\.csv$", full.names = TRUE)
  if (length(train_file) != 1)
    stop("No finished explorer results in ", results_dir, ".\n", explain_folder(describe_folder(results_dir)),
         "\nThe train CSV is written only after all TRAIN configurations finish.", call. = FALSE)
  test_file <- sub("_train\\.csv$", "_test.csv", train_file)
  if (!file.exists(test_file)) stop("Missing test file: ", test_file, call. = FALSE)
  out_dir <- file.path(results_dir, ANALYSIS_FOLDER)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)   # created if missing
  if (!dir.exists(out_dir)) stop("Cannot create ", out_dir, call. = FALSE)
  train_all <- read.csv(train_file, stringsAsFactors = FALSE)
  test_all <- read.csv(test_file, stringsAsFactors = FALSE)
  if (nzchar(EXCLUDE_LABELS)) {
    train_all <- train_all[!grepl(EXCLUDE_LABELS, train_all$label), ]
    test_all <- test_all[!grepl(EXCLUDE_LABELS, test_all$label), ]
  }
  family <- sub("^explore_(.*)_train\\.csv$", "\\1", basename(train_file))
  mode_label <- sub("^.*_", "", basename(normalizePath(results_dir)))

  param_cols <- setdiff(grep("^p_", names(train_all), value = TRUE), "p_effect")
  desc_cols <- grep("^d_", names(train_all), value = TRUE)
  varying <- function(v) length(unique(train_all[[v]])) > 1
  param_cols <- param_cols[sapply(param_cols, varying)]
  num_params <- param_cols[sapply(train_all[param_cols], is.numeric)]
  cat_params <- setdiff(param_cols, num_params)
  rate_cols <- grep("^p_rate", param_cols, value = TRUE)

  # Numeric parameters are logged when they span more than a factor of 8, then standardised
  # with statistics learned on the TRAIN set only.
  spec <- lapply(setNames(num_params, num_params), function(v) {
    x <- train_all[[v]]; use_log <- all(x > 0) && max(x) / min(x) > 8
    y <- if (use_log) log(x) else x
    list(log = use_log, mean = mean(y), sd = sd(y))
  })
  model_frame <- function(d) {
    out <- data.frame(rej = round(d$rejection_rate * d$successful), n_ok = d$successful)
    out$miss <- out$n_ok - out$rej
    for (v in num_params) {
      x <- d[[v]]; if (spec[[v]]$log) x <- log(x)
      out[[v]] <- (x - spec[[v]]$mean) / spec[[v]]$sd
    }
    for (v in cat_params) out[[v]] <- factor(d[[v]], levels = sort(unique(train_all[[v]])))
    out
  }
  tree_frame <- function(d, cols = c(param_cols, desc_cols)) {
    x <- d[, cols, drop = FALSE]
    for (v in intersect(cat_params, names(x))) x[[v]] <- factor(x[[v]], levels = sort(unique(train_all[[v]])))
    x
  }
  fml <- function(terms) as.formula(paste("cbind(rej, miss) ~", paste(terms, collapse = " + ")))
  auc <- function(score, y) {
    y <- as.logical(y); n1 <- sum(y); n0 <- sum(!y)
    if (n1 == 0 || n0 == 0) return(NA_real_)
    r <- rank(score); (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  importance <- function(fit, scope = NULL) {
    d1 <- suppressWarnings(if (is.null(scope)) drop1(fit, test = "F") else drop1(fit, scope = scope, test = "F"))
    terms <- rownames(d1)[rownames(d1) != "<none>"]
    keep <- rownames(d1) != "<none>"
    data.frame(term = terms, importance_pct = 100 * (d1$Deviance[keep] - fit$deviance) / fit$null.deviance,
               p_value = d1[["Pr(>F)"]][keep], stringsAsFactors = FALSE)
  }
  eval_condition <- function(cond, data) {
    m <- regmatches(cond, regexec("^([A-Za-z_0-9.]+)\\s*(>=|<|=)\\s*(.*)$", cond))[[1]]
    if (length(m) == 0) return(rep(NA, nrow(data)))
    x <- data[[m[2]]]; val <- trimws(m[4])
    switch(m[3], ">=" = as.numeric(x) >= as.numeric(val), "<" = as.numeric(x) < as.numeric(val),
           "=" = as.character(x) %in% strsplit(val, ",")[[1]])
  }

  labels <- sort(unique(train_all$label))
  drivers <- list(); interactions <- list(); rules <- list(); rules_obs <- list(); validation <- list()
  calibration <- list()
  for (lab in labels) {
    for (truth in c("null", "alt")) {
      pick <- function(d) d[d$label == lab & d$successful >= 20 & !is.na(d$rejection_rate) &
                             (if (truth == "null") d$p_effect == 0 else d$p_effect != 0), ]
      tr <- pick(train_all); te <- pick(test_all)
      if (nrow(tr) < MIN_CONFIGS || nrow(te) < 5) {
        if (truth == "null") validation[[lab]] <- data.frame(label = lab, n_train = nrow(tr),
          n_test = nrow(te), note = "too few configurations", stringsAsFactors = FALSE)
        next
      }
      mf_tr <- model_frame(tr); mf_te <- model_frame(te)
      all_terms <- c(num_params, cat_params)
      all_terms <- all_terms[vapply(all_terms, function(v) length(unique(mf_tr[[v]])) > 1, logical(1))]   # drop terms that do not vary in this subset
      fit <- tryCatch(suppressWarnings(glm(fml(all_terms), family = quasibinomial, data = mf_tr)),
                      error = function(e) NULL)
      if (is.null(fit)) next
      imp <- importance(fit)
      drivers[[paste(lab, truth)]] <- cbind(label = lab, truth = truth, imp)
      top <- imp$term[order(-imp$importance_pct)][seq_len(min(TOP_TERMS, nrow(imp)))]
      fit2 <- fit
      if (length(top) >= 2) {
        fit2_try <- tryCatch(suppressWarnings(glm(as.formula(paste("cbind(rej, miss) ~",
          paste(all_terms, collapse = " + "), "+ (", paste(top, collapse = " + "), ")^2")),
          family = quasibinomial, data = mf_tr)), error = function(e) NULL)
        if (!is.null(fit2_try)) {
          fit2 <- fit2_try
          labs <- attr(terms(fit2), "term.labels"); ints <- labs[grepl(":", labs)]
          if (length(ints)) {
            ii <- tryCatch(importance(fit2, ints), error = function(e) NULL)
            if (!is.null(ii)) interactions[[paste(lab, truth)]] <- cbind(label = lab, truth = truth,
              ii[order(-ii$importance_pct), ][seq_len(min(10, nrow(ii))), ])
          }
        }
      }
      pred <- tryCatch(suppressWarnings(predict(fit2, newdata = mf_te, type = "response")),
                       error = function(e) rep(NA_real_, nrow(te)))
      mae <- mean(abs(pred - te$rejection_rate), na.rm = TRUE)
      mae_const <- mean(abs(mean(tr$rejection_rate) - te$rejection_rate))
      if (truth == "alt") {
        validation[[paste(lab, "alt")]] <- data.frame(label = lab, truth = "alt", n_train = nrow(tr),
          n_test = nrow(te), mae_model = mae, mae_constant = mae_const,
          skill = 1 - mae / mae_const, stringsAsFactors = FALSE)
        next
      }
      calibration[[lab]] <- data.frame(label = lab, predicted = pred, observed = te$rejection_rate)
      harm_tr <- tr$rejection_rate > HARM_LEVEL; harm_te <- te$rejection_rate > HARM_LEVEL
      # With only a handful of harmful configurations an AUC or a tree is noise, so both are skipped.
      reliable <- sum(harm_tr) >= MIN_HARMFUL_TRAIN && sum(harm_te) >= MIN_HARMFUL_TEST
      auc_rate <- NA_real_
      if (reliable && length(rate_cols) && any(harm_tr) && any(!harm_tr)) {
        rd <- function(d) { x <- as.data.frame(lapply(d[rate_cols], identity))
          cbind(x, as.data.frame(lapply(x, function(v) v^2), col.names = paste0(names(x), "_sq"))) }
        rate_fit <- tryCatch(suppressWarnings(glm(harm ~ ., family = binomial,
                   data = cbind(rd(tr), harm = harm_tr))), error = function(e) NULL)
        if (!is.null(rate_fit)) auc_rate <- auc(suppressWarnings(predict(rate_fit,
                   newdata = rd(te), type = "response")), harm_te)
      }
      # Fits a small tree on the chosen columns; returns its AUC on the TEST set and its rules.
      run_tree <- function(cols) {
        tx_tr <- tree_frame(tr, cols); tx_te <- tree_frame(te, cols)
        tx_tr$harm <- factor(harm_tr, levels = c(FALSE, TRUE))
        tree <- rpart::rpart(harm ~ ., data = tx_tr, weights = tr$successful, method = "class",
          control = rpart::rpart.control(minbucket = max(10, round(0.04 * nrow(tx_tr))),
                                         cp = 0.005, maxdepth = 4, xval = 5))
        frm <- tree$frame; leaf <- frm$var == "<leaf>"
        # A tree with a single leaf found no rule and predicts a constant, so it has no meaningful AUC.
        out <- list(auc = if (sum(leaf) > 1) auc(predict(tree, tx_te, type = "prob")[, "TRUE"], harm_te) else NA_real_,
                    rules = NULL)
        if (sum(leaf) > 1) {
          nodes <- as.integer(rownames(frm))[leaf]
          paths <- rpart::path.rpart(tree, nodes = nodes, print.it = FALSE)
          prob_harm <- frm$yval2[leaf, 5]
          rows <- lapply(seq_along(nodes), function(i) {
            conds <- paths[[i]][-1]
            in_tr <- Reduce(`&`, lapply(conds, eval_condition, data = tx_tr[names(tx_tr) != "harm"]))
            in_te <- Reduce(`&`, lapply(conds, eval_condition, data = tx_te))
            type <- if (prob_harm[i] >= 0.5) "HARMFUL region" else if (prob_harm[i] <= 0.05) "SAFE region" else "mixed region"
            data.frame(label = lab, type = type, rule = paste(conds, collapse = "  AND  "),
              n_train = sum(in_tr), harm_share_train = mean(harm_tr[in_tr]),
              n_test = sum(in_te), harm_share_test = if (sum(in_te)) mean(harm_te[in_te]) else NA_real_,
              stringsAsFactors = FALSE)
          })
          out$rules <- do.call(rbind, rows)
        }
        out
      }
      auc_tree <- NA_real_; auc_tree_obs <- NA_real_
      if (reliable && any(harm_tr) && any(!harm_tr)) {
        full_cols <- c(param_cols, desc_cols)
        res <- run_tree(full_cols); auc_tree <- res$auc
        if (!is.null(res$rules)) rules[[lab]] <- res$rules
        obs_cols <- setdiff(full_cols, UNOBSERVABLE)
        if (length(obs_cols) < length(full_cols) && length(obs_cols) >= 2) {
          res_obs <- run_tree(obs_cols); auc_tree_obs <- res_obs$auc
          if (!is.null(res_obs$rules)) rules_obs[[lab]] <- res_obs$rules
        }
      }
      auc_glm <- if (reliable && any(harm_te) && any(!harm_te)) auc(pred, harm_te) else NA_real_
      validation[[lab]] <- data.frame(label = lab, truth = "null", n_train = nrow(tr), n_test = nrow(te),
        harm_share_train = mean(harm_tr), harm_share_test = mean(harm_te),
        auc_rate_only = auc_rate, auc_model = auc_glm, auc_tree = auc_tree,
        auc_tree_observable = auc_tree_obs, mae_model = mae, mae_constant = mae_const, skill = 1 - mae / mae_const,
        n_harm_train = sum(harm_tr), n_harm_test = sum(harm_te),
        note = if (reliable) NA_character_ else "too few harmful configurations for AUC and rules",
        stringsAsFactors = FALSE)
    }
  }
  stack <- function(l) if (length(l)) { o <- do.call(rbind, l); rownames(o) <- NULL; o } else NULL
  drivers <- stack(drivers); interactions <- stack(interactions); rules <- stack(rules); rules_obs <- stack(rules_obs)
  validation <- if (length(validation)) do.call(rbind, lapply(validation, function(v) {
    for (m in c("truth", "harm_share_train", "harm_share_test", "auc_rate_only", "auc_model",
                "auc_tree", "auc_tree_observable", "mae_model", "mae_constant", "skill", "n_harm_train", "n_harm_test", "note"))
      if (!m %in% names(v)) v[[m]] <- NA
    v[, c("label", "truth", "n_train", "n_test", "harm_share_train", "harm_share_test", "n_harm_train",
          "n_harm_test", "auc_rate_only", "auc_model", "auc_tree", "auc_tree_observable", "mae_model", "mae_constant",
          "skill", "note")]
  })) else NULL
  if (is.null(validation) || all(!is.na(validation$note) & validation$note == "too few configurations")) {
    cat(sprintf("Family: %s (%s): NOTHING TO ANALYSE. Each test needs at least %d null configurations in the train set, and this run has fewer.\n  A pilot run is only a quick check of the script; run the explorer with MODE <- \"medium\" (a few minutes) or \"full\" and analyse again.\n",
                family, mode_label, MIN_CONFIGS))
    return(invisible(NULL))
  }
  rownames(validation) <- NULL
  write_table <- function(x, name) if (!is.null(x)) write.csv(x, file.path(out_dir, name), row.names = FALSE)
  write_table(drivers, "analysis_drivers.csv"); write_table(interactions, "analysis_interactions.csv")
  write_table(rules, "analysis_rules.csv"); write_table(rules_obs, "analysis_rules_observable.csv"); write_table(validation, "analysis_validation.csv")

  # ---- images ---------------------------------------------------------------------------
  open_device <- function(file) {
    args <- list(filename = file, width = FIG_WIDTH_IN, height = FIG_HEIGHT_IN, units = "in", res = FIG_DPI)
    if (FIG_FORMAT == "tiff") args$compression <- "lzw"
    if (FIG_FORMAT == "jpeg") args$quality <- 95
    if (isTRUE(capabilities("cairo"))) args$type <- "cairo"
    do.call(get(FIG_FORMAT), args)
  }
  grid <- function(k) c(ceiling(k / 3), min(k, 3))
  auc_cols <- c("auc_rate_only", "auc_model", "auc_tree", "auc_tree_observable")
  null_val <- validation[!is.na(validation$truth) & validation$truth == "null" &
                           rowSums(!is.na(validation[, auc_cols])) > 0, ]
  drv <- if (!is.null(drivers)) drivers[drivers$truth == "null", ] else NULL
  if (!is.null(drv) && nrow(drv)) {
    labs <- unique(drv$label); g <- grid(length(labs))
    open_device(file.path(out_dir, sprintf("analysis_drivers.%s", FIG_FORMAT)))
    par(mfrow = g, mar = c(4, 9, 3, 1), oma = c(0, 0, 3, 0))
    for (lab in labs) {
      d <- drv[drv$label == lab, ]; d <- d[order(-d$importance_pct), ][seq_len(min(8, nrow(d))), ]
      barplot(rev(d$importance_pct), names.arg = rev(sub("^p_", "", d$term)), horiz = TRUE, las = 1,
              cex.names = 0.7, col = "#0072B2", main = lab, cex.main = 0.9,
              xlab = "Share of null deviance explained (%)")
    }
    mtext(sprintf("What drives Type I error: %s (%s)", family, mode_label), outer = TRUE, line = 1)
    dev.off()
  }
  if (nrow(null_val)) {
    open_device(file.path(out_dir, sprintf("analysis_validation_auc.%s", FIG_FORMAT)))
    m <- t(as.matrix(null_val[, auc_cols])); missing_bar <- is.na(m); m[missing_bar] <- 0
    longest <- max(strwidth(null_val$label, units = "inches", cex = 0.7))
    par(mar = c(longest / par("csi") + 2, 4, 3, 1))          # room for the longest label
    mids <- barplot(m, beside = TRUE, names.arg = null_val$label, las = 2, cex.names = 0.7, ylim = c(0, 1.12),
            col = c("#999999", "#0072B2", "#D55E00", "#009E73"), ylab = "AUC on fresh test configurations",
            main = sprintf("Can each view predict which designs are harmful? %s (%s)", family, mode_label))
    abline(h = 0.5, lty = 3)
    if (any(missing_bar)) text(mids[missing_bar], 0.02, "no rule found", srt = 90, adj = c(0, 0.5), cex = 0.6, col = "gray30")
    legend("top", c("Missing rate only", "Fitted model", "Tree (all parameters)", "Tree (observable only)"), horiz = TRUE,
           fill = c("#999999", "#0072B2", "#D55E00", "#009E73"), bty = "n", cex = 0.75)
    dev.off()
  }
  cal <- stack(calibration)
  if (!is.null(cal)) {
    labs <- unique(cal$label); g <- grid(length(labs))
    open_device(file.path(out_dir, sprintf("analysis_calibration.%s", FIG_FORMAT)))
    par(mfrow = g, mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
    for (lab in labs) {
      d <- cal[cal$label == lab, ]
      plot(d$predicted, d$observed, pch = 16, cex = 0.6, col = "#0072B2", xlim = c(0, 1), ylim = c(0, 1),
           xlab = "Predicted Type I error", ylab = "Observed Type I error", main = lab, cex.main = 0.9)
      abline(0, 1, lty = 3); abline(h = 0.05, lty = 3, col = "gray60")
    }
    mtext(sprintf("Predictions against fresh test configurations: %s (%s)", family, mode_label),
          outer = TRUE, line = 1)
    dev.off()
  }
  cat("Family:", family, "| null configurations analysed per label (train/test):",
      paste(sprintf("%s %d/%d", validation$label[!is.na(validation$n_test) & validation$truth %in% "null"],
                    validation$n_train[!is.na(validation$n_test) & validation$truth %in% "null"],
                    validation$n_test[!is.na(validation$n_test) & validation$truth %in% "null"]), collapse = "; "), "\n")
  cat("Wrote tables and images to", normalizePath(out_dir), "\n")
  invisible(list(drivers = drivers, interactions = interactions, rules = rules, validation = validation))
}
run_analysis <- function() {
  if (!FIG_FORMAT %in% c("png", "tiff", "jpeg"))
    stop("FIG_FORMAT must be png, tiff or jpeg", call. = FALSE)
  targets <- if (nzchar(RESULTS_DIR)) RESULTS_DIR else find_runs(OUTPUT_ROOT)
  cat("Analysing:", paste(basename(targets), collapse = ", "), "\n\n")
  out <- lapply(targets, analyze_explorer); names(out) <- basename(targets)
  invisible(out)
}
# The whole workflow runs inside run_analysis(), so a problem stops the run cleanly even if you
# paste line by line. Each analysis folder is created if it is missing.
analysis_results <- run_analysis()
