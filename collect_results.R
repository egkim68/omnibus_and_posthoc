# COLLECT THE RESULTS INTO ONE UPLOAD PACKAGE
# Base R only. Copies every image and table of the runs in the chosen mode (from each run's "plots",
# "analysis" and "paper_figures" folders) into ONE flat folder, with the run and test in every file
# name so that nothing collides, writes RUN_INFO.txt into it (toolkit version, mode, families, time, the
# runs collected, and key numbers), and makes ONE UNIQUELY NAMED zip:
#   <OUTPUT_ROOT>/UPLOAD_<mode>_<date>_<time>.zip
# Older UPLOAD zips of the same mode are moved to <OUTPUT_ROOT>/old_uploads, so the newest file in the
# main folder is always the right one. On Windows the folder opens with the new zip selected.
# Run by the pipeline as its last step; it can also be run alone with source(). The mode comes from the
# environment variable SIM_MODE (the pipeline sets it), otherwise COLLECT_MODE below ("all" = every run).
# The flat folder _ALL_<mode> is reused and emptied each time (Windows often blocks deleting and
# recreating a folder that Explorer has open). If the zip cannot be made (R needs a zip tool), the
# folder is still complete: right-click it and choose Send to > Compressed (zipped) folder.

COLLECT_MODE <- "medium"
IS_WINDOWS <- .Platform$OS.type == "windows"
OUTPUT_ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "N:/paper-simulation-explorer")

# Makes a real zip of the files in `dest`. Tries R's own zip() (needs a zip program), and on Windows the
# tools built into Windows 10 and later: PowerShell's Compress-Archive, then tar.exe. Every attempt is
# checked to be a real zip file (it must start with "PK"); returns the name of the method, or FALSE.
make_zip <- function(dest, zipfile) {
  is_zip <- function() file.exists(zipfile) && isTRUE(file.size(zipfile) > 0) &&
    identical(readBin(zipfile, "raw", 2), as.raw(c(0x50, 0x4b)))
  reset <- function() if (file.exists(zipfile)) unlink(zipfile)
  status <- suppressWarnings(tryCatch(
    utils::zip(zipfile, files = list.files(dest, full.names = TRUE), flags = "-jq9X"), error = function(e) 127L))
  if (identical(as.integer(status), 0L) && is_zip()) return("R's zip program")
  reset()
  if (IS_WINDOWS) {
    d <- normalizePath(dest, winslash = "\\"); z <- normalizePath(zipfile, winslash = "\\", mustWork = FALSE)
    cmd <- sprintf("Compress-Archive -Path '%s\\*' -DestinationPath '%s' -Force", d, z)
    suppressWarnings(tryCatch(system2("powershell", c("-NoProfile", "-NonInteractive", "-Command", shQuote(cmd, type = "cmd")),
                                      stdout = FALSE, stderr = FALSE), error = function(e) NULL))
    if (is_zip()) return("Windows PowerShell")
    reset()
    suppressWarnings(tryCatch(system2("tar", c("-a", "-c", "-f", shQuote(z, type = "cmd"), "-C", shQuote(d, type = "cmd"), "."),
                                      stdout = FALSE, stderr = FALSE), error = function(e) NULL))
    if (is_zip()) return("Windows tar")
    reset()
  }
  FALSE
}

# Key numbers of a run, written into RUN_INFO.txt so a result package can be checked at a glance.
key_numbers <- function(run_dir) {
  family <- sub("_(pilot|medium|full)$", "", basename(run_dir)); out <- character(0)
  tr <- file.path(run_dir, sprintf("explore_%s_train.csv", family)); te <- file.path(run_dir, sprintf("explore_%s_test.csv", family))
  if (family == "posthoc" && file.exists(tr) && file.exists(te)) {
    d <- rbind(read.csv(tr), read.csv(te)); d <- d[d$successful >= 20 & d$p_effect == 0 & grepl("^FWER", d$label), ]
    if (nrow(d)) { m <- tapply(d$rejection_rate, d$label, mean)
      out <- c(sprintf("post hoc: mean familywise error over %d designs with no true difference", length(unique(paste(d$set, d$config_id)))), sprintf("  %-42s %.3f", names(m), m)) }
  }
  if (family == "anova") {
    sel <- file.path(run_dir, "analysis", "selector_evaluation.csv")
    if (file.exists(sel)) { s <- read.csv(sel, stringsAsFactors = FALSE)
      role <- if ("role" %in% names(s)) s$role else s$kind
      out <- c("anova on held-out designs: mean Type I error / share harmful / power   [role]",
               sprintf("  %-32s %.3f  %.3f  %.3f   [%s]", s$label, s$test_typeI_mean, s$test_share_harmful, s$test_power, role)) }
  }
  conf <- file.path(run_dir, "analysis", sprintf("confirmation_%s.csv", if (family == "anova") "omnibus" else "posthoc"))
  if (file.exists(conf)) { h <- read.csv(conf, stringsAsFactors = FALSE)
    out <- c(out, sprintf("pre-stated hypotheses (%s): %s", family, h$stage[1]),
             sprintf("  %-4s %-9s %s | %s", sub(" .*$", "", h$hypothesis), h$verdict, substr(h$part, 1, 40), h$value))
  }
  out
}

collect_results <- function() {
  mode <- Sys.getenv("SIM_MODE", unset = COLLECT_MODE)
  if (!mode %in% c("pilot", "medium", "full", "all")) stop('The mode must be "pilot", "medium", "full" or "all"', call. = FALSE)
  if (!dir.exists(OUTPUT_ROOT)) stop("Output folder not found: ", OUTPUT_ROOT, ". Check that the drive exists, or set OUTPUT_ROOT.", call. = FALSE)
  dirs <- list.dirs(OUTPUT_ROOT, recursive = FALSE)
  pattern <- if (mode == "all") "_(pilot|medium|full)$" else paste0("_", mode, "$")
  runs <- dirs[grepl(pattern, basename(dirs)) & !grepl("^_ALL_", basename(dirs))]
  if (!length(runs)) stop("No ", mode, " runs found under ", OUTPUT_ROOT, ". Run the explorers (and the analysis and plots) first.", call. = FALSE)
  dest <- file.path(OUTPUT_ROOT, paste0("_ALL_", mode))
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  if (dir.exists(dest)) {
    old <- list.files(dest, full.names = TRUE, recursive = TRUE)
    if (length(old)) {
      removed <- suppressWarnings(file.remove(old))
      if (!all(removed)) cat(sprintf("Note: %d old file(s) in %s could not be removed (open in another program?).\n", sum(!removed), dest))
    }
  } else {
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(dest)) {
      dest <- paste0(dest, "_", stamp)
      dir.create(dest, recursive = TRUE, showWarnings = FALSE)
      if (!dir.exists(dest)) stop("Cannot create a collection folder under ", OUTPUT_ROOT, ". Check that the drive is connected and writable.", call. = FALSE)
      cat("The usual folder name could not be created, so this run uses ", dest, "\n", sep = "")
    }
  }
  copied <- 0L; failed <- 0L; per_run <- integer(0)
  for (run in runs) {
    n <- 0L
    for (sub in c("plots", "analysis", "paper_figures")) {
      from_dir <- file.path(run, sub)
      if (!dir.exists(from_dir)) next
      for (r in list.files(from_dir, pattern = "\\.(png|csv|tiff|jpeg)$", recursive = TRUE)) {
        to <- file.path(dest, paste0(basename(run), "__", gsub("/", "__", r, fixed = TRUE)))
        if (isTRUE(file.copy(file.path(from_dir, r), to, overwrite = TRUE))) n <- n + 1L else failed <- failed + 1L
      }
    }
    per_run[basename(run)] <- n; copied <- copied + n
  }
  cat(sprintf("Collected %d files from %d run(s) into %s\n", copied, length(runs), normalizePath(dest)))
  for (nm in names(per_run)) cat(sprintf("  %-24s %d files%s\n", nm, per_run[[nm]], if (per_run[[nm]] == 0L) "   <- nothing yet: run the analysis and plots for this run" else ""))
  if (failed) cat(sprintf("Note: %d file(s) could not be copied (open in another program?).\n", failed))

  # RUN_INFO.txt: what this package is, so it can be verified without opening anything else
  script_dir <- if (exists("SCRIPT_DIR")) SCRIPT_DIR else file.path(OUTPUT_ROOT, "R")
  version <- if (file.exists(file.path(script_dir, "VERSION.txt"))) readLines(file.path(script_dir, "VERSION.txt"), warn = FALSE)[1] else "unknown"
  fingerprint <- if (file.exists(file.path(script_dir, "MANIFEST.csv"))) unname(tools::md5sum(file.path(script_dir, "MANIFEST.csv"))) else "no manifest"
  info <- c("RUN INFO", "========",
            sprintf("Toolkit version:      %s (manifest fingerprint %s)", version, fingerprint),
            sprintf("Created:              %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
            sprintf("Mode:                 %s", mode),
            sprintf("Families run now:     %s", Sys.getenv("SIM_FAMILIES", unset = "(collector run alone)")),
            sprintf("R version:            %s", R.version.string), sprintf("Platform:             %s", R.version$platform),
            "Runs collected:", sprintf("  %-24s %d files", names(per_run), per_run), "",
            "Key numbers:")
  kn <- unlist(lapply(runs, key_numbers)); info <- c(info, if (length(kn)) kn else "  (none for these runs)")
  writeLines(info, file.path(dest, "RUN_INFO.txt"))

  # older UPLOAD zips of this mode are moved aside, so the newest file is always the right one
  old_up <- list.files(OUTPUT_ROOT, pattern = sprintf("^UPLOAD_%s_.*\\.zip$", if (mode == "all") "all" else mode), full.names = TRUE)
  if (length(old_up)) {
    hold <- file.path(OUTPUT_ROOT, "old_uploads"); dir.create(hold, showWarnings = FALSE)
    for (f in old_up) suppressWarnings(file.rename(f, file.path(hold, basename(f))))
  }
  zipfile <- file.path(OUTPUT_ROOT, sprintf("UPLOAD_%s_%s.zip", mode, stamp))
  method <- make_zip(dest, zipfile)
  if (!isFALSE(method)) {
    cat("\n=====================================================================\n")
    cat(sprintf("UPLOAD THIS FILE:  %s\n(%.1f MB, %d files, made with %s)\n", normalizePath(zipfile), file.size(zipfile) / 1024^2, copied + 1L, method))
    cat("=====================================================================\n")
    if (IS_WINDOWS && interactive())
      try(system2("explorer.exe", sprintf("/select,\"%s\"", gsub("/", "\\\\", normalizePath(zipfile, winslash = "/"))), wait = FALSE), silent = TRUE)
  } else {
    cat("\nCould not make a zip automatically (no zip tool worked).\n",
        "Right-click the folder ", normalizePath(dest), " and choose Send to > Compressed (zipped) folder, then upload that.\n", sep = "")
  }
  invisible(dest)
}
# The whole workflow runs inside collect_results(), so a problem stops the run cleanly.
collected_folder <- collect_results()
