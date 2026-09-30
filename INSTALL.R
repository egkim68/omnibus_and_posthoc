# INSTALL THE TOOLKIT (run once per release; see README.txt for the one line that starts it)
# Base R only. What it does, in order, and it stops with a plain message if any step fails:
#   1  checks that the downloaded toolkit is complete and undamaged (every file against MANIFEST.csv)
#   2  checks that the drive and folder exist
#   3  BACKS UP whatever is in the script folder to <root>/R_old/<date>_<time> (nothing is deleted)
#   4  installs every file and verifies each installed file against the manifest
#   5  moves stale result zips from earlier workflows (_ALL_*.zip, "files (N).zip") to <root>/old_uploads
#   6  runs the setup check (CHECK_setup.R) and tells you the next command
TARGET_DIR <- Sys.getenv("SIM_TARGET_DIR", unset = "N:/paper-simulation-explorer/R")     # where the scripts are installed
INSTALL_HERE <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")), error = function(e) NA_character_)

install_toolkit <- function(here, target) {
  if (is.na(here)) stop("Start the installer with the one line in README.txt (it uses source()).", call. = FALSE)
  manifest <- file.path(here, "MANIFEST.csv")
  if (!file.exists(manifest)) stop("MANIFEST.csv is not next to INSTALL.R. The toolkit download is incomplete; download it again.", call. = FALSE)
  m <- read.csv(manifest, stringsAsFactors = FALSE)
  src <- unname(tools::md5sum(file.path(here, m$file)))
  if (any(is.na(src) | src != m$md5))
    stop("The downloaded toolkit is damaged or incomplete (", paste(m$file[is.na(src) | src != m$md5], collapse = ", "), "). Download it again.", call. = FALSE)
  version <- readLines(file.path(here, "VERSION.txt"), warn = FALSE)[1]
  root <- dirname(target)
  if (!dir.exists(root)) stop("The folder ", root, " was not found. Check that the drive (N:) is connected, or change TARGET_DIR at the top of INSTALL.R.", call. = FALSE)
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(target)) stop("Cannot create ", target, call. = FALSE)
  cat(sprintf("Installing toolkit version %s into %s\n", version, target))
  # 3: back up whatever is there
  present <- list.files(target, full.names = TRUE); present <- present[!file.info(present)$isdir]
  if (length(present)) {
    backup <- file.path(root, "R_old", stamp); dir.create(backup, recursive = TRUE, showWarnings = FALSE)
    copied <- file.copy(present, backup, overwrite = TRUE)
    same <- unname(tools::md5sum(present)) == unname(tools::md5sum(file.path(backup, basename(present))))
    if (!all(copied) || any(is.na(same) | !same))
      stop("The backup to ", backup, " could not be verified, so nothing was removed or installed. ",
           "Check the free space and permissions of that folder and run the installer again.", call. = FALSE)
    file.remove(present)
    cat(sprintf("  backed up %d existing file(s) to %s\n", length(present), backup))
  }
  # 4: install and verify
  files <- c(m$file, "MANIFEST.csv")                 # the manifest itself is installed too: the run guard needs it
  ok <- file.copy(file.path(here, files), file.path(target, files), overwrite = TRUE)
  installed <- unname(tools::md5sum(file.path(target, m$file)))
  manifest_ok <- identical(unname(tools::md5sum(file.path(target, "MANIFEST.csv"))), unname(tools::md5sum(manifest)))
  if (!all(ok) || !manifest_ok || any(is.na(installed) | installed != m$md5))
    stop("Installation failed for: ", paste(c(m$file[is.na(installed) | installed != m$md5], if (!manifest_ok) "MANIFEST.csv", files[!ok]), collapse = ", "),
         ". Is a file open in another program? Close it and run the installer again.", call. = FALSE)
  cat(sprintf("  installed and verified %d files (plus the manifest)\n", nrow(m)))
  # 5: stale result zips from earlier workflows
  stale <- list.files(root, pattern = "^_ALL_.*\\.zip$|^files \\([0-9]+\\)\\.zip$", full.names = TRUE)
  if (length(stale)) {
    hold <- file.path(root, "old_uploads"); dir.create(hold, showWarnings = FALSE)
    moved <- 0L; for (f in stale) if (suppressWarnings(file.rename(f, file.path(hold, basename(f))))) moved <- moved + 1L
    cat(sprintf("  moved %d old zip file(s) to %s so they cannot be uploaded by mistake\n", moved, hold))
  }
  # 6: the setup check
  TOOLKIT_DIR <<- target
  source(file.path(target, "CHECK_setup.R"), local = FALSE)
  cat("\nNext: start the control panel with this one line, and follow its NEXT ACTION:\n")
  cat(sprintf("  source(\"%s/START.R\")\n", target))
  invisible(TRUE)
}
install_toolkit(INSTALL_HERE, TARGET_DIR)
