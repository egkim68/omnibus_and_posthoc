# CHECK THE SETUP: reports what works on this computer and what does not, with the fix for each.
# Run any time with  source("N:/paper-simulation-explorer/R/CHECK_setup.R")
# It changes nothing except a small temporary file it creates and removes in the output folder.

check_setup <- function() {
  root <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "N:/paper-simulation-explorer")
  dir <- if (exists("TOOLKIT_DIR")) TOOLKIT_DIR else tryCatch(dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")), error = function(e) file.path(root, "R"))
  rows <- list()
  add <- function(check, ok, detail = "", fix = "") rows[[length(rows) + 1L]] <<- data.frame(check = check, ok = ok, detail = detail, fix = fix, stringsAsFactors = FALSE)

  v <- getRversion(); add("R version is 4.0 or newer", v >= "4.0.0", as.character(v), "Install a current R from r-project.org")
  add("rpart package is available", requireNamespace("rpart", quietly = TRUE), "", "In R:  install.packages(\"rpart\")")
  add("output folder exists", dir.exists(root), root, "Connect the drive (N:) or create the folder")
  writable <- FALSE
  if (dir.exists(root)) { tmp <- file.path(root, ".write_test.tmp"); writable <- tryCatch({ writeLines("x", tmp); ok <- file.exists(tmp); unlink(tmp); ok }, error = function(e) FALSE, warning = function(w) FALSE) }
  add("output folder is writable", writable, root, "Check that the drive is connected and not read-only")
  # the installed scripts must be exactly the release
  mf <- file.path(dir, "MANIFEST.csv"); if (file.exists(mf)) {
    m <- read.csv(mf, stringsAsFactors = FALSE); actual <- unname(tools::md5sum(file.path(dir, m$file))); bad <- m$file[is.na(actual) | actual != m$md5]
    add("scripts match the installed release", !length(bad), if (length(bad)) paste(bad, collapse = ", ") else sprintf("%d files verified", nrow(m)), "Run the installer again")
  } else add("scripts match the installed release", FALSE, "MANIFEST.csv not found", "Run the installer again")
  # can a zip be made? test every method on a tiny folder
  zip_method <- "none"
  td <- tempfile("ziptest"); dir.create(td); writeLines("hello", file.path(td, "a.txt")); zf <- tempfile(fileext = ".zip")
  is_zip <- function() file.exists(zf) && isTRUE(file.size(zf) > 0) && identical(readBin(zf, "raw", 2), as.raw(c(0x50, 0x4b)))
  st <- suppressWarnings(tryCatch(utils::zip(zf, files = list.files(td, full.names = TRUE), flags = "-jq9X"), error = function(e) 127L))
  if (identical(as.integer(st), 0L) && is_zip()) zip_method <- "R's zip program"
  if (zip_method == "none" && .Platform$OS.type == "windows") {
    unlink(zf); d <- normalizePath(td, winslash = "\\"); z <- normalizePath(zf, winslash = "\\", mustWork = FALSE)
    suppressWarnings(tryCatch(system2("powershell", c("-NoProfile", "-NonInteractive", "-Command", shQuote(sprintf("Compress-Archive -Path '%s\\*' -DestinationPath '%s' -Force", d, z), type = "cmd")), stdout = FALSE, stderr = FALSE), error = function(e) NULL))
    if (is_zip()) zip_method <- "Windows PowerShell" else { unlink(zf)
      suppressWarnings(tryCatch(system2("tar", c("-a", "-c", "-f", shQuote(z, type = "cmd"), "-C", shQuote(d, type = "cmd"), "."), stdout = FALSE, stderr = FALSE), error = function(e) NULL))
      if (is_zip()) zip_method <- "Windows tar" }
  }
  unlink(c(td, zf), recursive = TRUE)
  add("a zip file can be made automatically", zip_method != "none", zip_method, "Not fatal: right-click the result folder > Send to > Compressed (zipped) folder")
  add("text encoding", TRUE, paste(Sys.getlocale("LC_CTYPE"), "(the scripts use plain ASCII, so any locale works)"))

  res <- do.call(rbind, rows)
  cat("\nSETUP CHECK\n===========\n")
  for (i in seq_len(nrow(res))) {
    cat(sprintf("  [%s] %s%s\n", if (res$ok[i]) "PASS" else "FAIL", res$check[i], if (nzchar(res$detail[i])) paste0("  (", res$detail[i], ")") else ""))
    if (!res$ok[i]) cat(sprintf("         fix: %s\n", res$fix[i]))
  }
  hard <- res$check[!res$ok & res$check != "a zip file can be made automatically"]
  cat(if (length(hard)) sprintf("\nNOT READY: %d problem(s) above need fixing.\n", length(hard)) else "\nREADY: everything needed is in place.\n")
  invisible(res)
}
check_setup_result <- check_setup()
if (exists("TOOLKIT_DIR")) suppressWarnings(rm(TOOLKIT_DIR, envir = globalenv()))
