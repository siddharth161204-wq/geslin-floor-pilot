# Independent recomputation

`recompute.py` is a second implementation of the headline quantities, written from the definitions in `predictions.md` without access to `pilot_a_floor.py` or its outputs. It verifies the input hashes, then recomputes the pair counts, the distribution of the pair difference at 90 and 85 % SOH, sigma, the resampled worst-of-8 and worst-of-96 shortfall (10 000 000 simulated modules each), the within-pair variance share, the worse-cell estimator dependence and the paper's Fig. 2c averages, and writes `recompute.json`. `crosscheck.py` executes the computational code cells of the authors' figure2.ipynb and compares their output with the re-implementation, checks the lifetime interpolation against `np.interp`, and repeats the worst-of-N medians with other seeds.

```
python3 -I verification/recompute.py dynamic_cycling_Nature_Energy_2024/data
python3 -I verification/crosscheck.py dynamic_cycling_Nature_Energy_2024/data
```

Agreement with `results/results.json`: identical pair counts (45 at 90 %, 43 at 85 %), mean, median, 95th percentile and maximum of the pair difference at both thresholds, sigma, within-pair variance share (0.029887) and estimator dependence (36 of 45 pairs); the nine Fig. 2c averages agree with the literal notebook run to 0.0; the worst-of-8 and worst-of-96 medians agree to Monte Carlo precision (4.116 % and 10.162 % here from 10 000 000 modules, 4.10 % and 10.16 % in the main script from 50 000).

The logs are `run_log.txt` and `crosscheck_log.txt`. The block `alt_B` in `run_log.txt` gives the pair statistics if the lifetime is instead the first crossing of the threshold on a curve filtered from the start; it is not part of the pre-registered analysis.
