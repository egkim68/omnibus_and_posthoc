# Omnibus and Post Hoc Tests Under Missing Data

Repository: https://github.com/egkim68/omnibus_and_posthoc

Simulation code for the manuscript *Omnibus and Post Hoc Tests Under Missing Data: A Simulation Study of Type I Error in Small, Unbalanced Groups* (Eungi Kim, Keimyung University).

The study simulates random one-way designs (3 to 40 groups) with missing data and compares five omnibus tests (F, Welch, Brown-Forsythe, Alexander-Govern, Kruskal-Wallis) and four post hoc procedures (Tukey-Kramer, Games-Howell, Welch-Holm, Dunn-Holm), plus four gated pipelines. An exploratory phase generated hypotheses; a confirmatory phase with a new declared seed tested them against prespecified criteria. All data, tables, and figures in the manuscript can be regenerated from this code and the declared seeds.

## Contents

| Folder | What it holds |
|---|---|
| `code/` | The simulation toolkit, release 1.8 (base R plus `rpart`). `MANIFEST.csv` fingerprints every file. |
| `docs/` | Prespecified hypotheses, the protocol, and the decision record. |

## Running the simulation

Requirements: R 4.0 or newer (the runs used R 4.3.3) and the `rpart` package, which ships with R.

1. Clone the repository (`git clone https://github.com/egkim68/omnibus_and_posthoc`).
2. Set the output folder, for example `Sys.setenv(SIM_OUTPUT_ROOT = "C:/sim")`, and copy the files of `code/` to `<output folder>/R`.
3. Run `source("<output folder>/R/START.R")` and follow its NEXT ACTION. It checks the setup, runs the exploratory ("medium") runs, and, after the rules are recorded as frozen, the confirmatory ("full") runs.
4. Each run ends with a zip whose `RUN_INFO.txt` lists the toolkit version and key numbers.

The omnibus confirmatory run was simulated with release 1.4 and re-reported with release 1.5 from the unchanged saved results; the post hoc confirmatory run used release 1.5. Release 1.6 corrected figures only. Releases 1.7 and 1.8 change only reporting, figures, installer safety, and comments, so they reproduce the same values. `code/reproduce_text_numbers.R` prints the numbers reported in the text that the pipeline does not print.

Seeds: exploratory 2033 (omnibus) and 2035 (post hoc); confirmatory 732033 and 732035. The exploratory runs take minutes; each confirmatory run takes about 1.5 hours on a desktop computer.

## Deviations

Two deviations are recorded in `docs/decision_record.md`: the reporting script of release 1.4 re-chose the omnibus threshold on the confirmatory training designs (corrected in release 1.5 without rerunning the simulation), and the scoring code for the criteria was written after the omnibus results were seen.

## Licence

Code: MIT licence (`LICENSE`). Documentation: CC BY 4.0.
