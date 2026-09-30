# Prespecified hypotheses and criteria

Stated before the confirmatory run (seeds 732033 and 732035). "Harmful" means an error rate above 0.10 in a design.

- **H1**: for at least 5 of 6 key procedures, the held-out AUC of the missing rate alone is below 0.70 and the AUC of the design-aware model is at least 0.15 higher.
- **H2**: the precision-only harmful share of Welch ANOVA exceeds that of Brown-Forsythe by at least 20 points, and Brown-Forsythe's bias-only share is at least 15% and at least 10 points above its precision-only share.
- **H3**: the omnibus rule's harmful share on all full designs is at most 15% and at least 10 points below both Welch ANOVA's and the F test's.
- **H4**: Welch-Holm's harmful share is at most 15% and at least 10 points below both Games-Howell's and Tukey-Kramer's.
- **H4b**: with outliers, Dunn-Holm's per-pair power exceeds Welch-Holm's by at least 0.05 with a harmful share of at most 15%; under equal means with skewed data and unequal spread, its harmful share exceeds Welch-Holm's.

Frozen rule: Brown-Forsythe when the smallest observed group has fewer than 5 observations, Welch ANOVA otherwise.

How each criterion is scored (implemented in `code/check_hypotheses.R`, release 1.5; overall H1 verdict added in 1.7 and made to require all 13 procedure rows in 1.8):

- H1 per procedure on held-out designs (`analysis_validation.csv`). In the confirmatory run, every procedure with an estimable AUC met the criterion and one was not estimable, so H1 holds for any choice of six key procedures.
- H2 on null designs, training and held-out pooled; precision channel = smallest expected group below 4; bias channel = outcome-driven missingness.
- H3 on all null designs, training and held-out pooled.
- H4 on null designs, pooled.
- H4b (a): designs with outliers and a true difference (per-pair power), and null designs with outliers (harmful share); (b): null designs with skewed data and SD ratio above 1.

The original protocol, including its list of the six key procedures for H1, is `PROTOCOL_v1.md` in this folder.
