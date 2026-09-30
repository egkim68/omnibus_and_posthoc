# Decision record

| Date | Decision or event | Reason and evidence |
|---|---|---|
| 2026-09-23 | Rules frozen after the exploratory phase: omnibus threshold 5 (Brown-Forsythe below, Welch otherwise); criteria H1 to H4b fixed. | Threshold chosen on exploratory training designs as the most powerful rule with at most 10% harmful designs. |
| 2026-09-30 | Omnibus confirmatory run completed with toolkit 1.4. | 1,500 training and 500 held-out designs, 200 replications, seed 732033. |
| 2026-09-30 | **Deviation 1.** The 1.4 reporting script re-chose the threshold on the confirmatory training designs (it selected 8) and reported that rule as the main result. | This contradicted the freeze. Toolkit 1.5 reports the frozen rule; the re-chosen threshold is kept as a sensitivity row. The reporting was rerun from the saved results; the simulation was not repeated and no simulated value changed. |
| 2026-09-30 | **Deviation 2.** The scoring code for the criteria (`check_hypotheses.R`) was written into toolkit 1.5 after the omnibus results had been seen. | It follows the prespecified wording, including "all full designs" for H3. The subsets that define H4b were fixed at this point, before the post hoc run. |
| 2026-09-30 | Post hoc confirmatory run completed with toolkit 1.5 (seed 732035). | Verdicts: H1, H2, H4, H4b met; H3 not met (11.2% harmful; 9.6 points below the F test, 10 required). |
| 2026-09-30 | Toolkit 1.6: figure corrections only (Figure 6 axis label, figure annotations and legends). | No results changed. |
| 2026-09-30 | Toolkit 1.7: reporting scripts follow the requested mode; installer verifies its backup; comments on target parameters. | No simulated value changed. |
| 2026-09-30 | Toolkit 1.8: H1 overall verdict requires all 13 procedure rows and is saved to a file; H4b (a) compares power on designs where both estimates are finite and reports that number; `reproduce_text_numbers.R` added for numbers not printed by the pipeline; Figure 5A shows only the rule family of the paper, and Figures 5B and 6B relabelled. | Review of the code and manuscript. No simulated value changed; the H4b (a) means are unchanged because a design lacks finite power for every procedure at once. |
