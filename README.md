# The cell-to-cell lifetime floor in the Geslin et al. dynamic-cycling dataset

How far apart do nominally identical lithium-ion cells drift in lifetime under one protocol at one temperature, by chance alone? Geslin et al. (Nature Energy, 2025) cycled 92 commercial SiOx-graphite/NCA cells under 47 protocols at 35 °C, each protocol on two cells, and report that the average lifetime difference between duplicates is below 5 %. A floor is a tail, not an average. This repository measures the distribution of the duplicate difference and converts it into a selection-only floor: how far below the module mean the worst of N identical cells falls by chance. A thermal explanation of cell-to-cell lifetime divergence in a module has to exceed this floor.

Author: Siddharth Satte. Archived at https://doi.org/10.5281/zenodo.23205413.

## Pre-registration

The predictions in [predictions.md](predictions.md) were committed before the analysis was run, in commit `409eba9326fbfc0fbcc5b81aba4c044b1fd54a21`, and pushed to this repository before the script was run on the data; the tag and pre-release `pre-registration` mark that commit. The analysis script in that commit is the one that produced the results. The first commit of this repository holds the 7 October 2026 draft of the script; every change between it and the pre-registered version is listed below, and all were made before any result was computed.

Release v1.0.1 changes names only: the working label used while the study was being planned was removed, and the analysis script and its figure were renamed. The file predictions.md differs from the pre-registered one only where it used that label or named the script (the title, the opening paragraph and the second gate); no definition, prediction, band or rule changed, and the tag `pre-registration` holds the original files. Rerunning the renamed script reproduces every file in `results/` byte for byte.

## Results

The one-page note is [NOTE.md](NOTE.md) (also [NOTE.pdf](NOTE.pdf)); every number is in [results/results.json](results/results.json).

- 8.7 %: the 95th percentile of the lifetime difference between duplicate cells at 90 % SOH (percentile bootstrap over pairs, 5.7 to 13.9 %). The mean is 3.0 % and the median 1.7 %; 9 of 45 pairs differ by more than 5 %.
- 4.1 %: the median shortfall of the worst of eight nominally identical cells below the module mean lifetime, by selection alone, at 90 % SOH (2.8 to 5.6 %); the normal model gives 3.6 %.
- All nine scored predictions held; the tail was heavier than predicted. The scorecard is in the note.

### Gates, checks and reconciliation

- Gates on the processing passed: 45 pairs at 90 % and 43 at 85 % SOH; the pair differences equal those of the paper's notebook processing within 0.0002 for every pair except Synthetic_2c_Co2 at 85 % SOH, which contains cell_085. There the notebook's single pass keeps three points of cell_085 that a later diagnostic exceeds, and the pair difference moves by 0.10.
- The paper's Fig. 2c group averages, recomputed with its own notebook code, are all below 5 % (largest 4.89 %, C/2 at 85 % SOH; 3.21 % with the monotone filter), as the paper states. The ratio of the range of the duplicate differences to the protocol spread is at most 0.50 at 87.5 and 85 % SOH, as the paper states for SOH beyond 90 %, and 0.53 (C/5) and 0.59 (C/2) at 90 % SOH itself.
- The within-pair share of lifetime variance at 90 % SOH is 3.0 %, against 97.0 % between protocols.
- Found while reconciling, not pre-registered: the authors' file of lifetimes at 90 % SOH (`eol_metrics_soh0.9.pkl`) agrees with the lifetimes computed here within 1 % for 88 of 92 cells. The other four (cell_035, cell_043, cell_075, cell_081) sit at the first crossing of the threshold on the unfiltered curve, up to 95 EFC earlier. With that file the pair difference at 90 % SOH has mean 3.8 %, 95th percentile 10.8 % and maximum 23.1 %.
- Sensitivity: S1 (the notebook's processing) is identical at 90 % SOH and gives a mean of 2.95 % and a 95th percentile of 10.46 % at 85 % SOH; S2 (without the C/16 drive protocols) gives a 95th percentile at 90 % SOH of 8.70 % and a worst-of-8 median of 3.88 %; S3 (the authors' analysis set, 42 pairs) gives 8.75 % and 3.82 %.
- An independent implementation in [verification/](verification/) reproduces every deterministic number exactly and the simulated ones to Monte Carlo precision.

There were no deviations from the pre-registered plan in the processing or the scoring. The run used Python 3.13.16, numpy 2.5.3, pandas 3.0.5 and matplotlib 3.11.2 (recorded in `results/results.json`).

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
python3 lifetime_floor.py --data dynamic_cycling_Nature_Energy_2024/data --out results --upstream-commit 5b2f7f04d05072fe1f9bd8af664f23eadbde6317
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

The CSV files hold the per-cell estimators, the pair differences, the worse-cell agreement matrix and the per-pair ratio ranges; `lifetime_floor.png` is the figure.

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

## Citation

Satte, S. (2026). The cell-to-cell lifetime floor in the Geslin et al. dynamic-cycling dataset (v1.0.1). Zenodo. https://doi.org/10.5281/zenodo.23205413

The archived version is the GitHub release v1.0.1 of this repository. Please also cite the Geslin et al. paper and dataset listed below.

## Licence

Code: MIT (LICENSE). Results, figure and note: CC BY 4.0. Data: Geslin et al., CC BY 4.0, not included.

## References

Geslin, A., Xu, L., Ganapathi, D., Moy, K., Chueh, W. C. & Onori, S. Dynamic cycling enhances battery lifetime. Nature Energy 10, 172-180 (2025). https://doi.org/10.1038/s41560-024-01675-8

Geslin, A., Xu, L., Ganapathi, D., Moy, K., Chueh, W. & Onori, S. Dataset - Dynamic cycling enhances battery lifetime. Stanford Digital Repository (2024). https://doi.org/10.25740/td676xr4322
