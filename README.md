# Pilot A: the cell-to-cell lifetime floor in the Geslin et al. dynamic-cycling dataset

How far apart do nominally identical lithium-ion cells drift in lifetime under one protocol at one temperature, by chance alone? Geslin et al. (Nature Energy, 2025) cycled 92 commercial SiOx-graphite/NCA cells under 47 protocols at 35 °C, each protocol on two cells, and report that the average lifetime difference between duplicates is below 5 %. A floor is a tail, not an average. This repository measures the distribution of the duplicate difference and converts it into a selection-only floor: how far below the module mean the worst of N identical cells falls by chance. A thermal explanation of cell-to-cell lifetime divergence in a module has to exceed this floor.

Author: Siddharth Satte.

## Pre-registration

The predictions in [predictions.md](predictions.md) were committed before the analysis was run, in commit `PREREG_COMMIT` (tagged `pre-registration`). The analysis script in that commit is the one that produced the results. The first commit of this repository holds the 7 October 2026 draft of the script; every change between it and the pre-registered version is listed below, and all were made before any result was computed.

## Data

The data are not redistributed here. The script reads the processed per-cell files that the paper's authors released in their code repository, https://github.com/geslina/dynamic_cycling_Nature_Energy_2024, at commit `5b2f7f04d05072fe1f9bd8af664f23eadbde6317`: `data/metadata.pkl` (cell, protocol, protocol type and variant, nominal and measured average C-rate) and `data/diagnostic_features_all.pkl` (C/2 discharge capacity, normalised, and EFC at each diagnostic). The reconciliation also reads `data/eol_metrics_soh0.9.pkl`. Expected SHA-256:

```
metadata.pkl                 7579f48cf3767552154fba3ad9acb048c3166088c2a46fadb06b2a221b3a1e1f
diagnostic_features_all.pkl  906e54db95cf9a3f294b8fa0ba29d43b5850a7f1a228cdd9e37cd2bcdfe1bc79
eol_metrics_soh0.9.pkl       04a7427f563a1f6be5ec7b13193ad1c9812a9a0fad6e4fbba691378e30859f48
```

The script checks these hashes before it unpickles anything and stops on a mismatch. It then loads the pickles through a whitelist of the exact pandas and numpy constructors they use; any other global is refused, so a substituted pickle cannot run code. The full dataset is at the Stanford Digital Repository under a CC BY 4.0 licence (doi.org/10.25740/td676xr4322).

## How to run

```
git clone https://github.com/geslina/dynamic_cycling_Nature_Energy_2024
git -C dynamic_cycling_Nature_Energy_2024 checkout 5b2f7f04d05072fe1f9bd8af664f23eadbde6317
python3 -m pip install numpy pandas matplotlib
python3 pilot_a_floor.py --data dynamic_cycling_Nature_Energy_2024/data --out results --upstream-commit 5b2f7f04d05072fe1f9bd8af664f23eadbde6317
python3 tests/synthetic_check.py --data dynamic_cycling_Nature_Energy_2024/data
```

The last command runs the script on synthetic trajectories with known noise built on the real experimental design and prints twelve PASS lines; it reads only `metadata.pkl` from the real data. Random draws use a fixed seed (20261007).

## What the script reports

`results.json` holds every number.

- A1: the duplicate-pair relative lifetime difference at 95, 92.5, 90, 87.5 and 85 % SOH (mean, median, 75th, 90th and 95th percentiles, maximum), overall, by protocol type and by nominal C-rate.
- A2: the per-cell noise sigma implied by the pair differences, from their mean and median with bootstrap intervals and from their root mean square, and tail diagnostics against the values the same number of normal pairs would give.
- A3: the selection-only worst-of-N shortfall for N = 2 to 96, from a normal model and by resampling the pair deviations; percentile-bootstrap intervals over pairs for the headline numbers (indicative only for the 95th percentile, whose upper limit cannot exceed the largest observed pair).
- A4: which duplicate is the worse cell under twelve lifetime and fade-rate estimators, the agreement matrix, and the range across estimators of the pair degradation ratio.
- A5: the within-pair and between-protocol shares of lifetime variance.
- R: the reconciliation with the paper's Fig. 2c and with its file of lifetimes at 90 % SOH.
- S1 to S3: the sensitivity runs described below.

The CSV files hold the per-cell estimators, the pair differences, the worse-cell agreement matrix and the per-pair ratio ranges; `pilot_a_floor.png` is the figure.

## Processing

The C/2 discharge capacity of each cell is made monotonically decreasing by dropping every point that a later point exceeds, which is the rule the paper states in its Methods ("Non-monotonically decreasing capacity data points were ignored to ensure reliable EoL criteria"). The lifetime at each threshold is the EFC at which the capacity falls to that fraction of the first retained capacity, by linear interpolation; a cell that never reaches a threshold has no lifetime there.

The paper's notebook (figure2.ipynb) applies the rule in a single backward pass that compares each point with its original predecessor, and then interpolates the whole SOH grid in one vectorised call. In this dataset the single pass leaves two curves (cell_077 and cell_085) non-monotone, and on a non-monotone curve the interval that numpy's search lands in can depend on the neighbouring grid point. The notebook's processing is therefore reproduced exactly as sensitivity S1 and for the Fig. 2c reconciliation, and is not the primary.

Sensitivity S2 leaves out the two C/16 drive protocols (cells 089, 090, 093 and 094), which the paper excluded from its analyses because they did not reach the lower cut-off voltage. Sensitivity S3 is the authors' own analysis set: their figure3.ipynb also drops cell_045, so S3 additionally leaves out the Periodic C pair at C/2, whose two cells have the largest gap in measured average C-rate in the design (0.611 against 0.599 per hour). The scorecard uses the primary analysis only.

The two cells whose duplicate failed, cell_017 (CC, no rest, C/2) and cell_084 (synthetic profile 2c, C/5), are excluded from pair statistics by construction; the two shorted cells, 018 and 083, are absent from the released files.

## Changes from the 7 October draft, all made before any result was computed

1. Filter. The draft's filter is kept as the primary and its description corrected: it is the paper's stated rule carried to a monotone curve, not the notebook's exact code. The notebook's single-pass filter and its vectorised interpolation are added verbatim for sensitivity S1 and the reconciliation.
2. Shortfall. The worst-of-N shortfall is computed as defined, (mean minus minimum) divided by the mean of the simulated lifetimes, instead of mean minus minimum of the deviations; the difference is second order. Simulated modules per N rise from 20 000 to 50 000.
3. Ratio. The pair degradation ratio is also reported oriented on the cell that is worse at 90 % SOH, with lifetimes entering inverted, so that it falls below 1 where an estimator disagrees; this is the analogue of a ratio that crosses 1. The draft's maximum-over-minimum ratio, which is always at least 1, is kept. The worse-cell flip fraction is also reported within the five lifetimes, within the seven rates and fades, and among the ten estimators that compare the two cells over matched windows (the two whole-record estimators, average fade rate to the last diagnostic and slope over all diagnostics, use each cell's own record). An exact tie between the two cells names no worse cell for that estimator.
4. Estimators. The average fade rate to the last diagnostic is measured from the first retained diagnostic, whose EFC is not zero. Fade at 200 and 400 EFC is undefined, rather than clamped, if the retained record starts after that EFC. The first and last EFC of each retained record are written to the per-cell table.
5. A1 by group. The per-group 95th percentile is reported only where a group has at least ten pairs, a per-group maximum is added, and the C-rate grouping is renamed by_nominal_Crate with labels C/16, C/10, C/5 and C/2.
6. A2. The interval for sigma from the median and sigma from the root mean square (without an interval) are added; the draft already reported sigma from the mean and the median.
7. Added: tail diagnostics; percentile-bootstrap intervals over pairs for the headline numbers; the variance split on logarithmic lifetimes as a secondary measure; the reconciliation; sensitivity runs S1 to S3; input hashes checked before loading and package versions recorded; the exact-name whitelist loader; a check that EFC increases strictly in every cell; non-finite values written to JSON as null; the normal model drawn in the figure next to the resampled one.
8. Added `tests/synthetic_check.py`.

## Licence

Code: MIT (LICENSE). Results, figure and note: CC BY 4.0. Data: Geslin et al., CC BY 4.0, not included.

## References

Geslin, A., Xu, L., Ganapathi, D., Moy, K., Chueh, W. C. & Onori, S. Dynamic cycling enhances battery lifetime. Nature Energy 10, 172-180 (2025). https://doi.org/10.1038/s41560-024-01675-8

Geslin, A., Xu, L., Ganapathi, D., Moy, K., Chueh, W. & Onori, S. Dataset - Dynamic cycling enhances battery lifetime. Stanford Digital Repository (2024). https://doi.org/10.25740/td676xr4322
