# Predictions, written before the analysis ran

Siddharth Satte, 7 October 2026. This file was committed before the analysis script was run on the data for results. The commit that adds it is recorded in README.md and tagged `pre-registration`. The analysis script in that commit is the one that produces the results.

## What had been seen when these predictions were written

1. The paper (Geslin et al., Nature Energy 10, 172-180, 2025). Its Fig. 2c shows, for the 13, 9 and 6 protocols cycled at C/10, C/5 and C/2, the average lifetime difference between duplicate cells at 90, 87.5 and 85 % SOH, and the text states that this average "is consistently below 5%" and that the range "does not explain more than half of the variability beyond 90% SOH". No value was read off the figure.
2. The code of the paper's repository at commit 5b2f7f04d05072fe1f9bd8af664f23eadbde6317, read without its stored notebook outputs, including the list of cells its figure3.ipynb excludes (cell_045, cell_089, cell_090, cell_093, cell_094) and the cell names in its training and test lists.
3. The structure of the data files, checked without computing any lifetime, pair difference or fade rate: 92 cells, 47 protocols, 45 duplicate pairs, singletons cell_017 (CC_(no_storage)_Co2) and cell_084 (Synthetic_2c_Co5), no missing values, EFC strictly increasing in every cell, 10 to 32 diagnostics per cell, and the measured average C-rate of every cell. One processing check counted the cells for which two forms of the monotonicity filter keep different points: two cells, cell_077 and cell_085. It computed no lifetime.
4. Normal-theory constants and simulations of 45 pairs under normal and heavier-tailed noise, which use no data.
5. The analysis script was executed on the data three times on 7 October 2026 (once as the draft, twice as the reviewed version) only to confirm that it runs; every output was deleted unseen each time. The reviewed script was also run on synthetic trajectories with known noise built on the real design (tests/synthetic_check.py), and was reviewed against this plan and the paper's notebook by a reader who saw no result.

## Definitions

Lifetime of a cell at a threshold: the EFC at which its C/2 discharge capacity, made monotonically decreasing by dropping every point that a later point exceeds, first falls to the threshold times its first retained capacity, by linear interpolation. Pair difference: |EFC_a - EFC_b| divided by the pair mean, at the same threshold for both cells. Per-cell noise sigma: the mean pair difference times sqrt(pi)/2. Worst-of-N shortfall: (mean minus minimum) divided by the mean of the lifetimes of N simulated cells; "resampled" draws each cell's deviation from the 90 pair deviations (plus and minus delta/2, scaled by sqrt 2), "normal" draws it from a normal distribution with the fitted sigma. Unless stated, values are at 90 % SOH over the 45 pairs (43 pairs reach 85 %).

## Predictions

Each prediction is a value with a band. It is scored held if the measured value lies inside the band and failed otherwise.

P1. Mean pair difference at 90 % SOH: 3.5 %. Band 2.5 to 5.0 %.
The paper bounds its C-rate group averages below 5 % at three thresholds; a bound stated as 5 % suggests averages of roughly 3 to 4.5 %. The pooled 45 pairs also contain the rest, C/16 and off-rate protocols that Fig. 2c leaves out, which I do not expect to move the mean much.

P2. 95th percentile of the pair difference: 9.0 % at 90 % SOH (band 6.5 to 12.5 %) and 9.5 % at 85 % SOH (band 6.5 to 13.5 %).
With the P1 mean, 45 pairs of normal noise give a 95th percentile near 8.2 %; the heavier tail expected below raises it towards 9. At 85 % SOH I expect duplicates to have drifted slightly further apart.

P3. Per-cell lifetime noise sigma at 90 % SOH: 3.1 % (P1 times sqrt(pi)/2). Band 2.2 to 4.4 %.

P4. Median selection-only shortfall of the worst cell at 90 % SOH, resampled: 4.2 % for N = 8 (band 3.0 to 6.0 %) and 8.5 % for N = 96 (band 6.0 to 12.5 %). The normal model with the P3 sigma gives 4.2 % and 7.6 %.
Worst-of-8 is about 1.37 sigma and depends little on the tail. With 45 pairs the resampled worst-of-96 is set by the largest pair difference, which appears in about two thirds of simulated modules of 96, so this part is in effect a prediction that the largest pair difference at 90 % SOH lies between 8.5 and 17.7 %. Worst-of-96 is an extrapolation beyond what 45 pairs can resolve and will be reported as one. Because a single anomalous cell enters the resampling pool at 1/sqrt 2 of its deviation, the resampled worst-of-N at large N can understate such a cell's effect by up to that factor.

P5. Fraction of pairs whose worse cell is not the same under all of the twelve estimators that are defined for both cells (five lifetimes, average fade rate to the last diagnostic, slopes over the first 4, 6 and 8 retained diagnostics and over all, fade at 200 and 400 EFC; at least three must be defined): 70 %. Band 50 to 90 %.
When two cells differ by about 3 % in lifetime, early slopes and fade at 200 EFC reflect formation-stage behaviour and diagnostic noise that need not follow the later ranking, and one disagreement among twelve estimators is enough to count. Two of the twelve, the average fade rate to the last diagnostic and the slope over all diagnostics, compare each cell over its own record, whose length can differ between duplicates; the fraction among the other ten is reported alongside and is not scored.

P6. Within-pair share of lifetime variance at 90 % SOH (EFC, not logarithm): 1.5 %. Band 0.4 to 4.0 %.
With per-cell noise near 3 % and a between-protocol coefficient of variation of 20 to 35 % (C-rate and rest effects plus the dynamic-cycling gains of up to 38 % the paper reports), the share sigma squared over the total is 0.8 to 2.4 %.

T. The tail is heavier than normal.
The 45 pairs pool 47 protocols whose per-cell noise need not be equal, and a mixture of normals with different widths has a heavier tail than any one of them; a single anomalous cell (a fixture contact, a channel, or a manufacturing defect short of failure) also enters as one outlying pair, and so can a capacity recovery that falls next to a threshold crossing. Scoring rule: held if and only if the ratio of the 95th percentile to the mean pair difference at 90 % SOH (p95_over_mean) exceeds the median of that ratio for 45 pairs of normal noise as computed in the script (p95_over_mean_normal_median, about 2.34). Expected ratio 2.6. The ratio of the maximum to the mean and the ratio of the median-based to the mean-based sigma are reported with it and are not scored. I do not expect the ratio to exceed the normal 95th percentile (about 2.80), so a held result is weak evidence on its own; that remark is not scored either.

## Fixed before the run, with no prediction attached

1. Primary processing: the monotone filter. Sensitivity S1 reproduces the paper notebook's processing exactly (single-pass filter, one vectorised interpolation over its SOH grid, EFC rounded to 0.1). Sensitivity S2 leaves out the two C/16 drive protocols (cells 089, 090, 093 and 094), which the paper excluded from its analyses because they did not reach the lower cut-off voltage. Sensitivity S3 is the authors' analysis set: S2 plus the Periodic C pair at C/2, because their figure3.ipynb drops cell_045; this pair also has the largest gap in measured average C-rate in the design (0.611 against 0.599 per hour). The scorecard uses the primary analysis only.
2. Gates on the processing, checked before any result is interpreted: 45 pairs at 90 % and 43 at 85 % SOH, and the primary pair differences equal the notebook's (within 0.001) for every pair without cell_077 or cell_085. If either fails, the processing is reconciled first and the reconciliation is recorded in the README as a deviation from this plan.
3. Checks reported as found, which do not change the processing whatever they show: the paper's Fig. 2c group averages recomputed with its own notebook logic, and whether all are below 5 %; whether the within-pair variance share is smaller than the between-protocol share.
4. Headline numbers: the 95th percentile of the pair difference at 90 % SOH and the median resampled worst-of-8 shortfall at 90 % SOH, each with a percentile-bootstrap interval over pairs (indicative only for the 95th percentile, whose upper limit cannot exceed the largest observed pair).
5. The inputs are checked against the released files by SHA-256 before anything is unpickled.
6. The singletons cell_017 and cell_084 are named and excluded from pair statistics by construction.
7. The result is the floor that a thermal explanation must exceed; it is not a test of any thermal hypothesis.
