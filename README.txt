SIMULATION TOOLKIT - README (version in VERSION.txt)

THE ONE COMMAND YOU NEED (after the installation below)
  source("N:/paper-simulation-explorer/R/START.R")
  It shows a status board (what is done, what is not, what is out of date), tells you the NEXT ACTION,
  and offers a numbered menu. Type 1 and press Enter to do the next action. It remembers what you have
  run, uploaded and frozen, keeps a history of every action, and will not start the confirmatory (full)
  run before the rules are frozen unless you insist.
  While R shows the busy cursor (hourglass) it is reading files on your drive, not waiting for you: the
  panel prints "reading the status..." first and reports how long the read took. Type only when the
  menu asks for a number.

INSTALL (once per release)
  1. Download the toolkit zip. Any file name is fine.
  2. In R, paste this ONE block, then choose the downloaded zip in the window that opens:

     z <- file.choose(); d <- file.path(tempdir(), "sim_toolkit"); unlink(d, recursive = TRUE); unzip(z, exdir = d, junkpaths = TRUE)
     for (inner in list.files(d, pattern = "[.]zip$", full.names = TRUE)) unzip(inner, exdir = d, junkpaths = TRUE)
     source(file.path(d, "INSTALL.R"))

  It checks the download, backs up your old scripts (nothing is deleted), installs, verifies every
  file, moves stale zips aside, and runs the setup check. It ends with READY or NOT READY.

WHEN A RUN FINISHES
  The console prints  UPLOAD THIS FILE: <path>  and Windows Explorer opens with that file selected.
  Upload exactly that file. Its name contains the date and time, so an old file cannot be mistaken
  for it. Inside is RUN_INFO.txt with the version, mode and key numbers. Then choose
  "I have uploaded the newest package" in the panel so it knows.

WHAT IS NEW IN 1.8 (reporting and figures only; the simulated values are unchanged)
  - check_hypotheses.R: the H1 overall verdict is given only when all 13 procedure rows (5 omnibus tests,
    4 post hoc procedures, 4 pipelines) are present, and it is written to confirmation_H1_overall_<mode>.csv.
    H4b (a) compares Dunn-Holm and Welch-Holm on the designs where both power estimates are finite and
    reports that number next to the number of eligible designs.
  - reproduce_text_numbers.R (new): prints the numbers in the text that the pipeline does not, among them the
    cap on group rates, slices, Table 2, the pipelines, the common subset, and the switching rule.
    Run it on its own:  Rscript reproduce_text_numbers.R  (uses SIM_MODE, default medium).
  - analyze_explorer.R reads the harmful line from SIM_HARM_LEVEL (default 0.10), for the sensitivity AUCs.
  - Figure 5 (make_paper_figures.R): panel A shows only the rule family used in the paper; panel B bins
    read "<3", "3 to <5", "5 to <10", "10 or more". Figure 6 (make_figure5.R): horizontal procedure labels
    and a legend key in the procedure colors.

WHAT IS NEW IN 1.7 (reporting and safety only; the simulated values are unchanged)
  - The reporting scripts use the run of the mode the pipeline was started in (SIM_MODE), so a medium
    pipeline no longer reports from an existing full run. Without SIM_MODE: full > medium > pilot.
  - check_hypotheses.R gives an overall H1 verdict: for the key procedures in H1_KEY_PROCEDURES if set,
    otherwise for any choice of six.
  - INSTALL.R verifies the backup before it removes the old scripts.
  - Comments only in explore_anova.R and explore_posthoc.R: p_N and the MAR rate are targets, the SD ratio
    is a spread parameter, and how undefined pair power is handled. Draws and results are unchanged.

WHAT IS NEW IN 1.6 (figures only)
  - Figure 6 (make_figure5.R): the x axis is correctly labelled "smallest expected group after missingness";
    tick labels and n= labels no longer overlap.
  - Figures 3 to 5 (make_paper_figures.R): one annotation when both trees find no rule; Figure 4 axis reads
    "Unique share of null deviance (%)"; the frozen-rule legend colour matches the plot.

WHAT IS NEW IN 1.5 (reporting only; the simulation code is unchanged)
  - The ANOVA selector table, Figure 5 and RUN_INFO.txt report the FROZEN rule (BF<5 else Welch) as the
    tested rule. A threshold re-chosen on the confirmatory data is shown only as a sensitivity row.
  - check_hypotheses.R writes confirmation_omnibus.csv and confirmation_posthoc.csv: each pre-stated
    criterion (H1 to H4b), the numbers behind it, and MET / NOT MET.
  - Figures 1 and 5: legends no longer cover data or labels. Figure 3 marks a tree with no rule.
    Figure_2_data.csv and Figure_5B_data.csv hold the exact numbers behind those figures.
  - For an ANOVA full run made with 1.4, the panel's NEXT ACTION is "REBUILD the ANOVA full results".
    It redoes only the reporting from the saved files (about 2 minutes). Nothing is simulated again.

IF SOMETHING GOES WRONG
  Copy the whole R console (click in it, Ctrl+A, Ctrl+C) and paste it in the chat.
  The scripts refuse to run if any file is not the installed version; the message says which file
  and how to fix it (run the installer again).

FILES (you never edit any of them)
  START.R                      the control panel
  RUN_<family>_<mode>.R        the launchers the panel uses (each can also be run alone with source())
  RUN_anova_full_reanalyse.R   rebuilds the ANOVA full reporting from saved results (no simulation)
  CHECK_setup.R                the setup check
  INSTALL.R, MANIFEST.csv      the installer and the fingerprints of every file
  explore_*.R, analyze_explorer.R, plot_explorer.R, evaluate_selectors.R, make_paper_figures.R,
  check_hypotheses.R, collect_results.R, run_pipeline.R   the procedure itself
