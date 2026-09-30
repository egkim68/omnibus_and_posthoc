# START: the control panel. The ONE command you need:
#   source("N:/paper-simulation-explorer/R/START.R")
# It shows a status board (what is done, what is not, what is out of date), tells you the NEXT ACTION,
# and offers a numbered menu. Choose 1 to do the next action; everything else is a manual override.
# It remembers what you have run, uploaded and frozen in <root>/project_state.csv and writes every
# action to <root>/run_history.csv (time, action, OK or FAILED, minutes). Base R only.
#
# Order of the procedure for Paper A (omnibus tests and post hoc procedures):
#   check setup > ANOVA medium > post hoc medium > upload the package > (Claude reviews) > rules frozen >
#   ANOVA full > post hoc full > upload the package.  The full runs are the confirmatory run: they are
#   done once, after the freeze, and the panel asks before starting one and refuses before the freeze
#   unless you insist.
# Version 1.5 adds one step for a full ANOVA run made with version 1.4 or earlier: "rebuild the ANOVA full
# results" redoes only the reporting from the saved results (no simulation), so that the frozen rule,
# not a threshold re-chosen on the confirmatory data, is the reported rule.

START_DIR <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")), error = function(e) "N:/paper-simulation-explorer/R")
ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "N:/paper-simulation-explorer")
STATE_FILE <- file.path(ROOT, "project_state.csv")
HISTORY_FILE <- file.path(ROOT, "run_history.csv")
CHOICE_QUEUE <- if (nzchar(Sys.getenv("SIM_START_CHOICES"))) as.integer(strsplit(Sys.getenv("SIM_START_CHOICES"), ",")[[1]]) else integer(0)   # for tests only

panel <- function() {
  # ---- guard: the scripts must be the installed release ------------------------------------------------
  manifest <- file.path(START_DIR, "MANIFEST.csv")
  if (!file.exists(manifest)) stop("MANIFEST.csv is missing from ", START_DIR, ". Run the installer again (see README.txt).", call. = FALSE)
  m <- read.csv(manifest, stringsAsFactors = FALSE)
  actual <- unname(tools::md5sum(file.path(START_DIR, m$file))); bad <- m$file[is.na(actual) | actual != m$md5]
  if (length(bad) && !identical(Sys.getenv("SIM_SKIP_CHECK"), "1"))
    stop("These files are missing or are not the installed version: ", paste(bad, collapse = ", "), ". Run the installer again (see README.txt).", call. = FALSE)
  version <- readLines(file.path(START_DIR, "VERSION.txt"), warn = FALSE)[1]
  if (!dir.exists(ROOT)) stop("The folder ", ROOT, " was not found. Check that the drive (N:) is connected.", call. = FALSE)

  # ---- state and history --------------------------------------------------------------------------------
  state_read <- function() if (file.exists(STATE_FILE)) read.csv(STATE_FILE, stringsAsFactors = FALSE, colClasses = "character") else data.frame(key = character(), value = character(), time = character(), stringsAsFactors = FALSE)
  state_get <- function(key) { s <- state_read(); v <- s$value[s$key == key]; if (length(v)) v[length(v)] else NA_character_ }
  state_set <- function(key, value) { s <- state_read(); s <- s[s$key != key, , drop = FALSE]
    write.csv(rbind(s, data.frame(key = key, value = value, time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), stringsAsFactors = FALSE)), STATE_FILE, row.names = FALSE) }
  log_action <- function(action, result, minutes, note = "") {
    h <- if (file.exists(HISTORY_FILE)) read.csv(HISTORY_FILE, stringsAsFactors = FALSE) else NULL
    write.csv(rbind(h, data.frame(time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), action = action, result = result, minutes = minutes, toolkit = version, note = substr(note, 1, 200), stringsAsFactors = FALSE)), HISTORY_FILE, row.names = FALSE) }

  # ---- what exists on disk (read ONCE per redraw: a slow drive makes every extra read expensive) ------------------
  inspect_run <- function(family, mode) {
    d <- file.path(ROOT, paste0(family, "_", mode))
    if (!dir.exists(d)) return(list(explorer = FALSE, when = as.POSIXct(NA), analysed = FALSE, stale = FALSE, saved = 0L))
    train <- file.path(d, sprintf("explore_%s_train.csv", family)); test <- file.path(d, sprintf("explore_%s_test.csv", family))
    explorer <- file.exists(train) && file.exists(test); when <- if (explorer) file.mtime(test) else as.POSIXct(NA)
    val <- file.path(d, "analysis", "analysis_validation.csv"); has_val <- explorer && file.exists(val); vt <- if (has_val) file.mtime(val) else as.POSIXct(NA)
    # the checkpoint folder holds hundreds of files, so it is listed ONLY for a run that is unfinished
    saved <- if (!explorer && dir.exists(file.path(d, "checkpoints"))) length(list.files(file.path(d, "checkpoints"))) else 0L
    list(explorer = explorer, when = when, analysed = has_val && vt >= when - 120, stale = has_val && vt < when - 120, saved = saved)
  }
  snapshot <- function() {
    t0 <- Sys.time(); cat(" (reading the status of your folders; on a slow drive this can take a few seconds)\n"); flush.console()
    runs <- list(); for (f in c("anova", "posthoc")) for (m in c("medium", "full")) runs[[paste0(f, "_", m)]] <- inspect_run(f, m)
    up <- list.files(ROOT, pattern = "^UPLOAD_.*[.]zip$", full.names = TRUE)                  # one listing for both modes
    latest <- function(mode) { z <- up[grepl(sprintf("/UPLOAD_%s_", mode), gsub("\\\\", "/", up))]; if (length(z)) z[which.max(file.mtime(z))] else NULL }
    snap <- list(runs = runs, up = list(medium = latest("medium"), full = latest("full")), state = state_read())
    snap$up_time <- lapply(snap$up, function(z) if (is.null(z)) as.POSIXct(NA) else file.mtime(z))
    snap$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs")); snap
  }
  sget <- function(snap, key) { v <- snap$state$value[snap$state$key == key]; if (length(v)) v[length(v)] else NA_character_ }
  run_when <- function(snap, mode) { w <- c(snap$runs[[paste0("anova_", mode)]]$when, snap$runs[[paste0("posthoc_", mode)]]$when); w <- w[!is.na(w)]; if (length(w)) max(w) else as.POSIXct(NA) }
  package_state <- function(snap, mode) {
    z <- snap$up[[mode]]; if (is.null(z)) return(list(text = "none yet", ok = FALSE, sent = FALSE, file = NULL))
    rw <- run_when(snap, mode); fresh <- is.na(rw) || snap$up_time[[mode]] >= rw - 120
    sent <- identical(sget(snap, paste0("sent_", mode)), basename(z))
    list(text = sprintf("%s  [%s]", basename(z), if (!fresh) "OUT OF DATE: a run finished after it was made" else if (sent) "sent" else "NOT SENT YET"), ok = fresh, sent = sent, file = z)
  }

  # ---- the status board -----------------------------------------------------------------------------------
  board <- function(snap) {
    cell_run <- function(r) if (r$explorer) paste("done", format(r$when, "%m-%d %H:%M")) else if (r$saved > 0) sprintf("INTERRUPTED (%d saved)", r$saved) else "not started"
    cell_an <- function(r) if (!r$explorer) "-" else if (r$analysed) "done" else if (r$stale) "OUT OF DATE" else "not done"
    w <- 24
    cat("\n", strrep("=", 74), "\n", sep = "")
    cat(sprintf(" PAPER A  STATUS BOARD      toolkit %s      %s\n", version, format(Sys.time(), "%Y-%m-%d %H:%M")))
    cat(strrep("=", 74), "\n", sep = "")
    cat(sprintf(" %-26s %-*s %-*s\n", "", w, "ANOVA (omnibus)", w, "Post hoc"))
    for (mode in c("medium", "full")) {
      a <- snap$runs[[paste0("anova_", mode)]]; p <- snap$runs[[paste0("posthoc_", mode)]]
      cat(sprintf(" %-26s %-*s %-*s\n", sprintf("%s run", if (mode == "medium") "Medium (exploration)" else "FULL (confirmatory)"), w, cell_run(a), w, cell_run(p)))
      cat(sprintf(" %-26s %-*s %-*s\n", "  analysis and plots", w, cell_an(a), w, cell_an(p)))
      cat(sprintf(" %-26s %s\n", sprintf("  %s package", mode), package_state(snap, mode)$text))
      if (mode == "medium") cat(sprintf(" %-26s %s\n", "Rules frozen", if (identical(sget(snap, "frozen"), "yes")) paste("yes,", substr(snap$state$time[snap$state$key == "frozen"], 1, 16)) else "no"))
    }
    cat(strrep("-", 74), "\n", sep = "")
    cat(sprintf(" (status read in %.1f seconds%s)\n", snap$seconds, if (snap$seconds > 8) "; your drive is slow, this is normal" else ""))
  }

  # ---- what to do next -------------------------------------------------------------------------------------
  next_action <- function(snap) {
    if (!identical(sget(snap, "setup_ok"), version)) return(list(id = "check", text = "check the setup (about 10 seconds)"))
    for (spec in list(c("anova", "medium", "ANOVA (medium), about 10 minutes"), c("posthoc", "medium", "post hoc (medium), about 10 minutes"))) {
      r <- snap$runs[[paste0(spec[1], "_", spec[2])]]
      if (!r$explorer || !r$analysed) return(list(id = paste0(spec[1], "_", spec[2]), text = paste("run", spec[3])))
    }
    pm <- package_state(snap, "medium")
    if (!pm$ok) return(list(id = "package_medium", text = "make the upload package for the medium runs (about 20 seconds)"))
    if (!pm$sent) return(list(id = "sent", text = sprintf("upload  %s  to the chat, then choose 'I have uploaded it'", basename(pm$file))))
    if (!identical(sget(snap, "frozen"), "yes")) return(list(id = "wait", text = "WAIT for Claude's review. When you both agree the rules are frozen, choose 'Record that the rules are frozen'"))
    a_full <- snap$runs[["anova_full"]]
    if (a_full$explorer && a_full$analysed && !identical(sget(snap, "anova_full_rebuilt"), version))
      return(list(id = "anova_full_rebuild", text = "REBUILD the ANOVA full results from the saved files (about 2 minutes, no simulation)"))
    pf0 <- package_state(snap, "full")
    if (identical(sget(snap, "anova_full_rebuilt"), version) && !snap$runs[["posthoc_full"]]$explorer && pf0$ok && !pf0$sent)
      return(list(id = "sent", text = sprintf("upload  %s  (the rebuilt ANOVA results) to the chat, then choose 'I have uploaded it'", basename(pf0$file))))
    for (spec in list(c("anova", "full", "ANOVA FULL (confirmatory), about 1.5 hours"), c("posthoc", "full", "post hoc FULL (confirmatory), about 1.5 hours"))) {
      r <- snap$runs[[paste0(spec[1], "_", spec[2])]]
      if (!r$explorer || !r$analysed) return(list(id = paste0(spec[1], "_", spec[2]), text = paste("run", spec[3])))
    }
    pf <- package_state(snap, "full")
    if (!pf$ok) return(list(id = "package_full", text = "make the upload package for the full runs"))
    if (!pf$sent) return(list(id = "sent", text = sprintf("upload  %s  to the chat, then choose 'I have uploaded it'", basename(pf$file))))
    list(id = "done", text = "all runs are complete. Next: the manuscript figures and tables, with Claude")
  }

  # ---- actions -----------------------------------------------------------------------------------------------
  ask <- function(title, options) {
    if (interactive()) return(utils::menu(options, title = title))
    if (length(CHOICE_QUEUE)) { v <- CHOICE_QUEUE[1]; CHOICE_QUEUE <<- CHOICE_QUEUE[-1]; return(v) }
    0L
  }
  run_script <- function(script, label, before = NULL) {
    t0 <- Sys.time(); msg <- ""
    ok <- tryCatch({ if (!is.null(before)) before(); source(file.path(START_DIR, script), local = FALSE, echo = FALSE); TRUE }, error = function(e) { msg <<- conditionMessage(e); FALSE })
    mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
    log_action(label, if (ok) "OK" else "FAILED", mins, msg)
    if (ok) cat(sprintf("\nDone in %.1f minutes.\n", mins)) else cat("\nThat step FAILED. Copy the whole console (click in it, Ctrl+A, Ctrl+C) and paste it in the chat.\n")
    ok
  }
  do_action <- function(id, snap) {
    if (id == "check") { ok <- run_script("CHECK_setup.R", "setup check")
      if (ok && exists("check_setup_result")) { r <- get("check_setup_result"); if (all(r$ok | r$check == "a zip file can be made automatically")) state_set("setup_ok", version) else cat("The check found problems (see above). Fix them, then choose the check again.\n") }
      return(invisible()) }
    if (id %in% c("anova_medium", "posthoc_medium", "anova_full", "posthoc_full")) {
      if (grepl("_full$", id)) {
        if (!identical(sget(snap, "frozen"), "yes")) {
          if (ask("\nThe rules are NOT frozen yet. The full run is the confirmatory run: it should be done ONCE, after the freeze.", c("Go back", "Run it anyway (not recommended)")) != 2) { cat("Not started.\n"); return(invisible()) } }
        if (ask("\nA full run takes about 1.5 hours on this computer. Keep the computer awake. If it is interrupted, choose this again and it resumes.", c("Yes, start it now", "No, go back")) != 1) { cat("Not started.\n"); return(invisible()) } }
      run_script(sprintf("RUN_%s.R", id), sprintf("run %s", id)); return(invisible()) }
    if (id == "anova_full_rebuild") {
      if (ask("\nRebuild the ANOVA full tables, figures and hypothesis checks from the saved results? No simulation is run; the results themselves do not change.", c("Yes, rebuild now", "No, go back")) != 1) { cat("Not started.\n"); return(invisible()) }
      if (run_script("RUN_anova_full_reanalyse.R", "rebuild anova_full (reporting only)")) state_set("anova_full_rebuilt", version)
      return(invisible()) }
    if (id %in% c("package_medium", "package_full")) { mode <- sub("^package_", "", id); run_script("collect_results.R", paste("package", mode), before = function() Sys.setenv(SIM_MODE = mode)); return(invisible()) }
    if (id == "sent") {
      done <- character(0)
      for (mode in c("medium", "full")) { z <- snap$up[[mode]]; if (!is.null(z) && !identical(sget(snap, paste0("sent_", mode)), basename(z))) { state_set(paste0("sent_", mode), basename(z)); done <- c(done, basename(z)) } }
      if (length(done)) { cat("Recorded as uploaded:", paste(done, collapse = ", "), "\n"); log_action("marked as uploaded", "OK", 0, paste(done, collapse = ", ")) } else cat("There is no new package to record.\n")
      return(invisible()) }
    if (id == "freeze") {
      if (ask("\nRecord that the rules are frozen? Do this ONLY after Claude and you have agreed the medium results are understood and the decisions are final.", c("Yes, record it", "No, go back")) == 1) { state_set("frozen", "yes"); log_action("rules frozen", "OK", 0); cat("Recorded. The full (confirmatory) runs are now allowed.\n") }
      return(invisible()) }
    if (id == "wait") { cat("\nNothing to run now. Send the medium package to Claude and wait for the review.\n"); return(invisible()) }
    if (id == "done") { cat("\nAll runs are complete.\n"); return(invisible()) }
  }

  # ---- the loop ------------------------------------------------------------------------------------------------
  repeat {
    snap <- snapshot(); board(snap); na <- next_action(snap)
    cat(" NEXT ACTION: ", na$text, "\n", sep = ""); cat(strrep("=", 74), "\n", sep = "")
    options <- c(paste("DO THE NEXT ACTION:", na$text), "Run ANOVA (omnibus), medium", "Run post hoc, medium", "Run ANOVA (omnibus), FULL (confirmatory)", "Run post hoc, FULL (confirmatory)",
                 "Check the setup", "I have uploaded the newest package", "Record that the rules are frozen", "Make the upload package again (medium)",
                 "Rebuild the ANOVA full results from saved files (no simulation)")
    pick <- ask("\nWHAT WOULD YOU LIKE TO DO?  (enter a number; just press Enter to quit)", options)
    if (pick == 0) { cat("Goodbye. Run START again any time.\n"); break }
    id <- switch(pick, na$id, "anova_medium", "posthoc_medium", "anova_full", "posthoc_full", "check", "sent", "freeze", "package_medium", "anova_full_rebuild")
    do_action(id, snap)
  }
  invisible(TRUE)
}
panel()
