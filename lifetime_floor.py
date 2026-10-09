#!/usr/bin/env python3
"""
The cell-to-cell lifetime floor in the Geslin et al. (Nature Energy 2025) dataset.

What it computes, from the processed per-cell diagnostics released with the paper
(github.com/geslina/dynamic_cycling_Nature_Energy_2024, data/):

  A1  Duplicate-pair relative lifetime difference at five SOH thresholds,
      delta = |EFC_a - EFC_b| / mean(EFC_a, EFC_b), with its distribution (mean, median,
      75th, 90th, 95th percentiles, maximum), overall, by protocol type and by nominal C-rate.
  A2  The per-cell lifetime noise sigma implied by the pair differences, from the mean and the
      median of delta with bootstrap intervals over pairs and from its root mean square, and tail
      diagnostics against the values the same number of normal pairs would give.
  A3  The selection-only "worst of N" shortfall, (mean - minimum) / mean of the lifetimes of
      N nominally identical cells, by chance alone: median and 95th percentile for
      N = 2, 4, 8, 12, 24, 48, 96, parametric (normal) and nonparametric (resampled pair
      deviations, scaled by sqrt 2 because a pair mean removes half the variance).
  A4  Estimator sensitivity: which duplicate is the worse cell under twelve lifetime and
      fade-rate estimators; the fraction of pairs whose worse cell depends on the estimator;
      the agreement matrix; and the range across estimators of the pair degradation ratio,
      oriented on the cell that is worse at 90 % SOH (the analogue of a ratio that crosses 1).
  A5  Variance split at each threshold: within-pair versus between-protocol share.
  R   Reconciliation with the paper: its Fig. 2c cell-to-cell averages and ranges recomputed
      with the paper's own notebook logic, and the paper's lifetimes at 90 % SOH from
      eol_metrics_soh0.9.pkl.
  S   Sensitivity: S1 the lifetimes exactly as the paper's notebook computes them; S2 without
      the two C/16 drive protocols that the paper set aside; S3 the authors' analysis set, which
      also drops the Periodic C pair at C/2 (cell_045 is excluded in their figure3.ipynb).

Processing. The capacity trajectory is made monotonically decreasing by dropping every point
that a later point exceeds (the paper's rule, Methods: "Non-monotonically decreasing capacity
data points were ignored to ensure reliable EoL criteria"), then the EFC at each SOH, relative
to the first retained diagnostic, is found by linear interpolation; a cell that never reaches a
threshold is NaN there. The paper's notebook (figure2.ipynb, code cell 10) applies the rule in
a single backward pass that compares each point with its original predecessor, then interpolates
the whole SOH grid in one call; in this dataset that pass leaves two curves (cell_077, cell_085)
non-monotone, where the interpolation becomes unreliable, so the notebook's processing is run as
sensitivity S1 and for the reconciliation, not as the primary.

Inputs are read only; nothing in the Geslin repository is modified. Their SHA-256 is checked
against the released files before anything is unpickled, and the pickles are loaded through a
whitelist of the exact pandas and numpy constructors they use; any other global is refused.

Usage:
  python3 lifetime_floor.py --data /path/to/dynamic_cycling_Nature_Energy_2024/data --out results
"""
import argparse
import hashlib
import importlib
import json
import math
import os
import pickle
import platform

import numpy as np
import pandas as pd

SEED = 20261007
CAPA = "C/2 discharge CC Capacity [normalized]"
EFC = "EFCs"
THRESHOLDS = [0.95, 0.925, 0.90, 0.875, 0.85]
N_LIST = [2, 4, 8, 12, 24, 48, 96]
N_BOOT = 5000          # bootstrap replicates for the sigma intervals
N_MC = 50000           # simulated modules per N for the worst-of-N point estimates
N_BOOT_HEAD = 2000     # bootstrap replicates for the headline intervals
N_MC_HEAD = 4000       # simulated modules per N inside each headline replicate
N_NULL = 20000         # simulated normal datasets for the tail reference values
REFERENCE = "life_soh0.900"        # orients the pair degradation ratio in A4
LIFE_COLS = [f"life_soh{s:.3f}" for s in THRESHOLDS]
RATE_COLS = ["rate_avg_last", "rate_slope_first4", "rate_slope_first6", "rate_slope_first8",
             "rate_slope_all", "fade_at_200efc", "fade_at_400efc"]
# The two protocols whose duplicate failed (shorted cells 018 and 083 are already absent).
SINGLETON_PROTOCOLS = ("CC_(no_storage)_Co2", "Synthetic_2c_Co5")
# Geslin et al., Methods: the two C/16 drive profiles (four cells) did not reach the lower
# cut-off voltage and "were excluded from subsequent analyses but kept in the dataset".
C16_DRIVE_PROTOCOLS = ("Drive_City1_Co16", "Drive_City2_Co16")
# The authors' figure3.ipynb (code cell 15) also drops cell_045, which leaves its duplicate cell_046
# without a pair; this pair has the largest gap in measured average C-rate of the 45 (0.611 and 0.599).
AUTHORS_ANALYSIS_SET_EXCLUDED = C16_DRIVE_PROTOCOLS + ("Periodic_C_Co2",)
WHOLE_RECORD_RATES = ("rate_avg_last", "rate_slope_all")   # each cell over its own record
CRATE_LABEL = {"Co16": "C/16", "Co10": "C/10", "Co5": "C/5", "Co2": "C/2"}
Z975 = 1.959963984540054     # standard normal 97.5th percentile


# ---------------------------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------------------------
# The exact globals the released files use (first block), and those pandas 3 writes for the
# synthetic check (second block). Every one is a data constructor; nothing else is resolved.
ALLOWED_GLOBALS = {
    ("builtins", "slice"), ("numpy", "dtype"), ("numpy", "ndarray"),
    ("numpy.core.multiarray", "_reconstruct"), ("numpy.core.multiarray", "scalar"), ("numpy.core.numeric", "_frombuffer"),
    ("pandas._libs.internals", "_unpickle_block"), ("pandas.core.frame", "DataFrame"),
    ("pandas.core.indexes.base", "Index"), ("pandas.core.indexes.base", "_new_Index"),
    ("pandas.core.indexes.multi", "MultiIndex"), ("pandas.core.indexes.range", "RangeIndex"),
    ("pandas.core.internals.managers", "BlockManager"),
    ("numpy._core.multiarray", "_reconstruct"), ("numpy._core.multiarray", "scalar"), ("numpy._core.numeric", "_frombuffer"),
    ("pandas", "DataFrame"), ("pandas", "Index"), ("pandas", "MultiIndex"), ("pandas", "StringDtype"),
    ("pandas.arrays", "StringArray"), ("pandas._libs.arrays", "__pyx_unpickle_NDArrayBacked"),
}
# SHA-256 of the released files at commit 5b2f7f04d05072fe1f9bd8af664f23eadbde6317.
EXPECTED_SHA256 = {
    "metadata.pkl": "7579f48cf3767552154fba3ad9acb048c3166088c2a46fadb06b2a221b3a1e1f",
    "diagnostic_features_all.pkl": "906e54db95cf9a3f294b8fa0ba29d43b5850a7f1a228cdd9e37cd2bcdfe1bc79",
    "eol_metrics_soh0.9.pkl": "04a7427f563a1f6be5ec7b13193ad1c9812a9a0fad6e4fbba691378e30859f48",
}


class _Whitelist(pickle.Unpickler):
    """Resolve only the exact (module, name) pairs in ALLOWED_GLOBALS; refuse everything else."""

    def find_class(self, module, name):
        if (module, name) not in ALLOWED_GLOBALS:
            raise pickle.UnpicklingError(f"refused global {module}.{name} in a data pickle")
        if module.startswith("numpy.core."):      # written by numpy 1.x; numpy 2.x renamed it
            try:
                alt = "numpy._core." + module[len("numpy.core."):]
                importlib.import_module(alt)
                module = alt
            except ImportError:
                pass
        return super().find_class(module, name)


def read_pickle(path):
    with open(path, "rb") as f:
        return _Whitelist(f).load()


def verify_inputs(data_dir, allow_other):
    """Check the inputs are the released files before anything is unpickled."""
    report = {}
    for fname, expected in EXPECTED_SHA256.items():
        path = os.path.join(data_dir, fname)
        if not os.path.exists(path):
            report[fname] = None
            continue
        got = sha256(path)
        report[fname] = got
        if got != expected and not allow_other:
            raise SystemExit(f"{fname}: SHA-256 {got} does not match the released file ({expected}). "
                             "Use the Geslin et al. repository at commit 5b2f7f0, or pass --allow-unverified-inputs.")
    return report


def monotone_filter(cycle, capacity):
    """Primary. Walk back from the last diagnostic; when a point is higher than the one before it,
    drop the earlier point and re-check the new neighbour, so the retained curve never rises."""
    c = np.array(cycle, dtype=float).copy()
    q = np.array(capacity, dtype=float).copy()
    i = len(c) - 1
    while i > 0:
        if q[i] > q[i - 1]:
            c = np.delete(c, i - 1)
            q = np.delete(q, i - 1)
            i = min(i, len(c) - 1)
        else:
            i -= 1
    return c, q


def notebook_filter(cycle, capacity):
    """Verbatim logic of the paper's figure2.ipynb (code cell 10): one backward pass comparing each
    point with its ORIGINAL predecessor; the earlier point is dropped when the later one is higher."""
    cycle = np.array(cycle, dtype=float)
    capacity = np.array(capacity, dtype=float)
    fc, fq = cycle.copy(), capacity.copy()
    for i in range(len(cycle) - 1, 0, -1):
        if capacity[i] > capacity[i - 1]:
            fc = np.delete(fc, i - 1)
            fq = np.delete(fq, i - 1)
    return fc, fq


def load(data_dir, filt):
    meta = read_pickle(os.path.join(data_dir, "metadata.pkl"))
    dfa = read_pickle(os.path.join(data_dir, "diagnostic_features_all.pkl"))
    raw, traj = {}, {}
    for cell in meta["cell_name"]:
        sub = dfa.loc[cell]
        e = np.asarray(sub[EFC], dtype=float)
        q = np.asarray(sub[CAPA].values, dtype=float)
        if np.isnan(e).any() or np.isnan(q).any() or np.any(np.diff(e) <= 0):
            raise ValueError(f"{cell}: missing values or EFC not strictly increasing")
        raw[cell] = (e, q)
        traj[cell] = filt(e, q)
    return meta, raw, traj


# ---------------------------------------------------------------------------------------------
# Per-cell estimators
# ---------------------------------------------------------------------------------------------
def efc_at_soh(efc, capa, thresholds):
    """Linear interpolation of EFC at SOH * capa[0]; NaN if the threshold is never reached."""
    capa0 = capa[0]
    out = []
    for s in thresholds:
        target = s * capa0
        if capa.min() > target:
            out.append(np.nan)
            continue
        out.append(float(np.interp(target, np.flip(capa), np.flip(efc), left=np.nan, right=np.nan)))
    return np.array(out)


def estimators(e, q):
    """Twelve ways to say how fast a cell degrades. Lifetimes: lower = worse. Rates and fades:
    higher = worse; they are in percent of the first retained capacity, per EFC or at an EFC."""
    out = {}
    for s, v in zip(THRESHOLDS, efc_at_soh(e, q, THRESHOLDS)):
        out[f"life_soh{s:.3f}"] = v
    soh = q / q[0]
    span = e[-1] - e[0]
    out["rate_avg_last"] = 100.0 * (1.0 - soh[-1]) / span if span > 0 else np.nan
    for k in (4, 6, 8):
        if len(e) >= k:
            A = np.vstack([e[:k], np.ones(k)]).T
            m, _ = np.linalg.lstsq(A, soh[:k], rcond=None)[0]
            out[f"rate_slope_first{k}"] = -100.0 * m
        else:
            out[f"rate_slope_first{k}"] = np.nan
    A = np.vstack([e, np.ones(len(e))]).T
    m, _ = np.linalg.lstsq(A, soh, rcond=None)[0]
    out["rate_slope_all"] = -100.0 * m
    for efc_fixed in (200.0, 400.0):
        if e[0] <= efc_fixed <= e[-1]:
            out[f"fade_at_{int(efc_fixed)}efc"] = 100.0 * (1.0 - float(np.interp(efc_fixed, e, soh)))
        else:
            out[f"fade_at_{int(efc_fixed)}efc"] = np.nan
    out["efc_first"] = float(e[0])     # descriptive, not estimators: the span of the retained record
    out["efc_last"] = float(e[-1])
    return out


def build_estimators(meta, traj):
    cells = list(meta["cell_name"])
    est = pd.DataFrame({c: estimators(*traj[c]) for c in cells}).T
    est.index.name = "cell"
    return est.join(meta.set_index("cell_name")[["protocol_name", "protocol_type", "avg_Crate"]])


# ---------------------------------------------------------------------------------------------
# Pairs and pair differences
# ---------------------------------------------------------------------------------------------
def make_pairs(meta, exclude=()):
    groups = meta.groupby("protocol_name")["cell_name"].apply(list)
    if any(len(v) > 2 for v in groups):
        raise ValueError("a protocol has more than two cells; the duplicate design assumed here does not hold")
    unknown = [p for p in exclude if p not in groups.index]
    if unknown:
        raise ValueError(f"excluded protocols not in the metadata: {unknown}")
    pairs = [(p, v[0], v[1]) for p, v in groups.items() if len(v) == 2 and p not in exclude]
    singletons = {v[0]: p for p, v in groups.items() if len(v) == 1}
    return groups, pairs, singletons


def pair_table(est, pairs):
    rows = []
    for s in THRESHOLDS:
        col = f"life_soh{s:.3f}"
        for p, a, b in pairs:
            la, lb = float(est.loc[a, col]), float(est.loc[b, col])
            if np.isnan(la) or np.isnan(lb):
                continue
            rows.append({"threshold": s, "protocol": p, "cell_a": a, "cell_b": b, "life_a": la, "life_b": lb,
                         "delta_rel": abs(la - lb) / ((la + lb) / 2.0),
                         "protocol_type": est.loc[a, "protocol_type"],
                         "nominal_Crate": CRATE_LABEL.get(est.loc[a, "avg_Crate"], est.loc[a, "avg_Crate"])})
    return pd.DataFrame(rows)


def deltas(pt, s):
    return pt.loc[np.isclose(pt["threshold"], s), "delta_rel"].to_numpy(dtype=float)


def summary(d):
    d = np.asarray(d, dtype=float)
    return {"n_pairs": int(len(d)), "mean": float(d.mean()), "median": float(np.median(d)),
            "p75": float(np.percentile(d, 75)), "p90": float(np.percentile(d, 90)),
            "p95": float(np.percentile(d, 95)), "max": float(d.max())}


# ---------------------------------------------------------------------------------------------
# A1 to A5
# ---------------------------------------------------------------------------------------------
def a1_summary(pt):
    out = {}
    for s in THRESHOLDS:
        sub = pt[np.isclose(pt["threshold"], s)]
        if sub.empty:
            continue
        summ = summary(sub["delta_rel"])
        for key, name in (("protocol_type", "by_protocol_type"), ("nominal_Crate", "by_nominal_Crate")):
            groups = {}
            for g, gg in sub.groupby(key):
                dd = gg["delta_rel"].to_numpy(dtype=float)
                groups[str(g)] = {"n": int(len(dd)), "mean": float(dd.mean()), "median": float(np.median(dd)),
                                  "max": float(dd.max()),
                                  "p95": float(np.percentile(dd, 95)) if len(dd) >= 10 else None}
            summ[name] = groups
        out[f"soh{s:.3f}"] = summ
    return out


def sigma_estimates(d):
    # For normal per-cell noise: E|eps_a - eps_b| = 2 sigma / sqrt(pi); median = 0.6745 sqrt(2) sigma;
    # E[(eps_a - eps_b)^2] = 2 sigma^2.
    return {"sigma_from_mean": float(d.mean() * math.sqrt(math.pi) / 2.0),
            "sigma_from_median": float(np.median(d) / (0.6745 * math.sqrt(2.0))),
            "sigma_from_rms": float(math.sqrt(np.mean(d ** 2) / 2.0))}


def a2_sigma(pt, rng):
    out = {}
    for s in THRESHOLDS:
        d = deltas(pt, s)
        if len(d) < 5:
            continue
        res = sigma_estimates(d)
        bm, bmed = np.empty(N_BOOT), np.empty(N_BOOT)
        for i in range(N_BOOT):
            db = d[rng.integers(0, len(d), len(d))]
            bm[i] = db.mean() * math.sqrt(math.pi) / 2.0
            bmed[i] = np.median(db) / (0.6745 * math.sqrt(2.0))
        res["sigma_ci95"] = [float(np.percentile(bm, 2.5)), float(np.percentile(bm, 97.5))]
        res["sigma_from_median_ci95"] = [float(np.percentile(bmed, 2.5)), float(np.percentile(bmed, 97.5))]
        res["n_pairs"] = int(len(d))
        out[f"soh{s:.3f}"] = res
    return out


def tail_diagnostics(d, rng):
    """Observed tail ratios against the same ratios for len(d) pairs of normal noise."""
    n = len(d)
    z = np.abs(rng.standard_normal((N_NULL, n)) - rng.standard_normal((N_NULL, n)))
    r95 = np.percentile(z, 95, axis=1) / z.mean(axis=1)
    rmax = z.max(axis=1) / z.mean(axis=1)
    obs95 = float(np.percentile(d, 95) / d.mean())
    obsmax = float(d.max() / d.mean())
    sig = sigma_estimates(d)
    return {"n_pairs": int(n),
            "p95_over_mean": obs95,
            "p95_over_mean_normal_median": float(np.median(r95)),
            "p95_over_mean_normal_95th_percentile": float(np.percentile(r95, 95)),
            "p95_over_mean_fraction_of_normal_datasets_at_or_above": float(np.mean(r95 >= obs95)),
            "max_over_mean": obsmax,
            "max_over_mean_normal_median": float(np.median(rmax)),
            "max_over_mean_normal_95th_percentile": float(np.percentile(rmax, 95)),
            "max_over_mean_fraction_of_normal_datasets_at_or_above": float(np.mean(rmax >= obsmax)),
            "sigma_from_median_over_sigma_from_mean": sig["sigma_from_median"] / sig["sigma_from_mean"],
            "pairs_above_2.772_sigma_from_mean": int(np.sum(d > Z975 * math.sqrt(2.0) * sig["sigma_from_mean"])),
            "pairs_expected_above_under_normal": 0.05 * n}


def worst_of_n_shortfall(eps):
    """(mean - min) / mean of the lifetimes 1 + eps of each simulated module (one per row)."""
    life = 1.0 + eps
    m = life.mean(axis=1)
    return (m - life.min(axis=1)) / m


def deviation_pool(d):
    # Each cell's deviation from its pair mean is +-delta/2; the pair mean removes half the
    # variance, so scale by sqrt(2) to recover a single-cell deviation.
    return np.concatenate([d / 2.0, -d / 2.0]) * math.sqrt(2.0)


def a3_worst(pt, a2, rng):
    out = {}
    for s in (0.90, 0.85):
        key = f"soh{s:.3f}"
        if key not in a2:
            continue
        d = deltas(pt, s)
        sig = a2[key]["sigma_from_mean"]
        pool = deviation_pool(d)
        per_n = {}
        for n in N_LIST:
            sp = worst_of_n_shortfall(rng.normal(0.0, sig, size=(N_MC, n)))
            snp = worst_of_n_shortfall(rng.choice(pool, size=(N_MC, n), replace=True))
            per_n[str(n)] = {"parametric_median": float(np.median(sp)), "parametric_p95": float(np.percentile(sp, 95)),
                             "nonparametric_median": float(np.median(snp)), "nonparametric_p95": float(np.percentile(snp, 95)),
                             "probability_the_pool_extreme_is_drawn": float(1.0 - (1.0 - 1.0 / len(pool)) ** n)}
        out[key] = {"sigma_used": float(sig), "pool_size": int(len(pool)),
                    "largest_single_cell_deviation_in_pool": float(d.max() / math.sqrt(2.0)),
                    "worst_of_N_shortfall": per_n}
    return out


def headline_bootstrap(d90, d85, rng):
    """Percentile bootstrap over pairs (2.5 and 97.5 percentiles) for the headline numbers. The
    simulated modules use common random numbers (the pool is sorted, so a fixed uniform draw maps
    to a fixed quantile of each resampled pool), so the interval reflects the pairs, not the MC.
    The maximum and the resampled worst-of-96 are left out: a resample contains the sample maximum
    with probability 0.64, so their percentile intervals are not confidence intervals. For the 95th
    percentile, the upper limit cannot exceed the largest observed pair."""
    z8, z96 = rng.standard_normal((N_MC_HEAD, 8)), rng.standard_normal((N_MC_HEAD, 96))
    u8 = rng.random((N_MC_HEAD, 8))
    keys = ["mean_90", "p95_90", "sigma_90", "np_worst_of_8_median_90",
            "par_worst_of_8_median_90", "par_worst_of_96_median_90", "mean_85", "p95_85"]
    acc = {k: np.empty(N_BOOT_HEAD) for k in keys}
    for i in range(N_BOOT_HEAD):
        b = d90[rng.integers(0, len(d90), len(d90))]
        sig = b.mean() * math.sqrt(math.pi) / 2.0
        pool = np.sort(deviation_pool(b))
        acc["mean_90"][i] = b.mean()
        acc["p95_90"][i] = np.percentile(b, 95)
        acc["sigma_90"][i] = sig
        acc["np_worst_of_8_median_90"][i] = np.median(worst_of_n_shortfall(pool[(u8 * len(pool)).astype(int)]))
        acc["par_worst_of_8_median_90"][i] = np.median(worst_of_n_shortfall(sig * z8))
        acc["par_worst_of_96_median_90"][i] = np.median(worst_of_n_shortfall(sig * z96))
        c = d85[rng.integers(0, len(d85), len(d85))]
        acc["mean_85"][i] = c.mean()
        acc["p95_85"][i] = np.percentile(c, 95)
    return {k: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))] for k, v in acc.items()}


def a4_estimators(est, pairs):
    """Worse cell and pair ratio under each estimator. The oriented ratio puts the cell that is
    worse at 90 % SOH in the numerator of the degradation speed (lifetimes enter inverted), so it
    is at least 1 for the reference and falls below 1 where another estimator disagrees."""
    worse, rows = {}, []
    for p, a, b in pairs:
        w, unoriented, oriented = {}, {}, {}
        ra, rb = est.loc[a, REFERENCE], est.loc[b, REFERENCE]
        ref_worse = None if (np.isnan(ra) or np.isnan(rb)) else (a if ra < rb else b)
        other = None if ref_worse is None else (b if ref_worse == a else a)
        for c in LIFE_COLS:
            x, y = est.loc[a, c], est.loc[b, c]
            if np.isnan(x) or np.isnan(y) or x == y:     # an exact tie names no worse cell
                continue
            w[c] = a if x < y else b
            unoriented[c] = max(x, y) / min(x, y)
            if ref_worse is not None:
                oriented[c] = est.loc[other, c] / est.loc[ref_worse, c]
        for c in RATE_COLS:
            x, y = est.loc[a, c], est.loc[b, c]
            if np.isnan(x) or np.isnan(y) or min(x, y) <= 0 or x == y:
                continue
            w[c] = a if x > y else b
            unoriented[c] = max(x, y) / min(x, y)
            if ref_worse is not None:
                oriented[c] = est.loc[ref_worse, c] / est.loc[other, c]
        worse[p] = w
        life_ids = {w[c] for c in LIFE_COLS if c in w}
        rate_ids = {w[c] for c in RATE_COLS if c in w}
        matched_ids = {v for c, v in w.items() if c not in WHOLE_RECORD_RATES}
        rows.append({"protocol": p, "cell_a": a, "cell_b": b, "worse_at_soh90": ref_worse,
                     "efc_last_a": est.loc[a, "efc_last"], "efc_last_b": est.loc[b, "efc_last"],
                     "n_estimators": len(w),
                     "n_agreeing_with_soh90": int(sum(1 for v in w.values() if v == ref_worse)),
                     "worse_cell_depends_on_estimator": len(set(w.values())) > 1,
                     "worse_cell_depends_among_lifetimes": len(life_ids) > 1,
                     "worse_cell_depends_among_rates": len(rate_ids) > 1,
                     "worse_cell_depends_among_window_matched": len(matched_ids) > 1,
                     "ratio_min": min(unoriented.values()) if unoriented else np.nan,
                     "ratio_max": max(unoriented.values()) if unoriented else np.nan,
                     "oriented_ratio_min": min(oriented.values()) if oriented else np.nan,
                     "oriented_ratio_max": max(oriented.values()) if oriented else np.nan})
    rr = pd.DataFrame(rows)
    names = LIFE_COLS + RATE_COLS
    agree = pd.DataFrame(index=names, columns=names, dtype=float)
    for c1 in names:
        for c2 in names:
            both = [p for p in worse if c1 in worse[p] and c2 in worse[p]]
            agree.loc[c1, c2] = np.mean([worse[p][c1] == worse[p][c2] for p in both]) if both else np.nan
    ok = rr[rr["n_estimators"] >= 3]
    res = {
        "n_estimators_defined": len(names),
        "n_pairs_with_3plus_estimators": int(len(ok)),
        "fraction_of_pairs_whose_worse_cell_depends_on_estimator": float(ok["worse_cell_depends_on_estimator"].mean()),
        "fraction_depends_among_the_five_lifetimes": float(ok["worse_cell_depends_among_lifetimes"].mean()),
        "fraction_depends_among_the_seven_rates_and_fades": float(ok["worse_cell_depends_among_rates"].mean()),
        "fraction_depends_among_the_ten_window_matched_estimators": float(ok["worse_cell_depends_among_window_matched"].mean()),
        "note_whole_record_estimators": "rate_avg_last and rate_slope_all use each cell's own record; "
                                        "the ten-estimator fraction leaves them out",
        "median_relative_gap_in_record_end_efc": float(
            (abs(ok["efc_last_a"] - ok["efc_last_b"]) / ((ok["efc_last_a"] + ok["efc_last_b"]) / 2)).median()),
        "pair_ratio_range_median_min": float(ok["ratio_min"].median()),
        "pair_ratio_range_median_max": float(ok["ratio_max"].median()),
        "pair_ratio_range_overall": [float(ok["ratio_min"].min()), float(ok["ratio_max"].max())],
        "oriented_ratio_median_min": float(ok["oriented_ratio_min"].median()),
        "oriented_ratio_median_max": float(ok["oriented_ratio_max"].median()),
        "oriented_ratio_overall": [float(ok["oriented_ratio_min"].min()), float(ok["oriented_ratio_max"].max())],
        "fraction_of_pairs_whose_oriented_ratio_crosses_1": float((ok["oriented_ratio_min"] < 1.0).mean()),
        "agreement_life90_vs_life85": float(agree.loc["life_soh0.900", "life_soh0.850"]),
        "agreement_life90_vs_slope_first6": float(agree.loc["life_soh0.900", "rate_slope_first6"]),
        "agreement_life90_vs_each_estimator": {c: float(agree.loc["life_soh0.900", c]) for c in names},
    }
    return res, agree, rr


def a5_variance(est):
    out = {}
    for s in THRESHOLDS:
        col = f"life_soh{s:.3f}"
        vals = est[[col, "protocol_name"]].dropna()
        twos = [v.to_numpy(dtype=float) for _, v in vals.groupby("protocol_name")[col] if len(v) == 2]
        if len(twos) < 5:
            continue
        res = {"n_pairs": len(twos)}
        for label, f in (("", lambda x: x), ("_log", np.log)):
            within = np.mean([np.var(f(v), ddof=1) for v in twos])
            total = np.var(f(np.concatenate(twos)), ddof=1)
            res[f"within_pair_variance_share{label}"] = float(within / total)
            res[f"between_protocol_variance_share{label}"] = float(1.0 - within / total)
        out[f"soh{s:.3f}"] = res
    return out


# ---------------------------------------------------------------------------------------------
# Reconciliation with the paper
# ---------------------------------------------------------------------------------------------
def notebook_lifetimes(meta, raw):
    """Every cell's EFC on the SOH grid np.linspace(1, 0.85, 31), exactly as figure2.ipynb (code cell
    10) computes it: single-pass filter, one vectorised np.interp call over the grid (on a non-monotone
    curve the interval numpy's search lands in can depend on the previous grid point), EFC rounded to
    0.1. Columns are the grid rounded to three decimals, as in the notebook."""
    sohs = np.linspace(1, 0.85, 31)
    rows = {}
    for cell in meta["cell_name"]:
        e, q = notebook_filter(*raw[cell])
        eol = np.round(np.interp(sohs * q[0], np.flip(q), np.flip(e), left=np.nan, right=np.nan), 1)
        eol[0] = 0
        rows[cell] = eol
    return pd.DataFrame.from_dict(rows, orient="index", columns=np.round(sohs, 3))


def paper_fig2c(meta, all_sohs, pt):
    """The paper's Fig. 2c cell-to-cell averages ('hatched areas') and ranges ('whiskers'), and the
    'maximum EFC difference with constant current' bars, recomputed with the logic of figure2.ipynb
    (code cells 8, 10, 11, 17) from notebook_lifetimes: buckets of 13, 9 and 6 protocols."""
    prots = [p for p in meta["protocol_name"].unique() if p not in SINGLETON_PROTOCOLS]
    intrinsic = {}
    for p in prots:
        two = meta.loc[meta["protocol_name"] == p, "cell_name"].tolist()
        m = all_sohs.loc[two].mean()
        m.iloc[0] = 1
        intrinsic[p] = (all_sohs.loc[two].max() - all_sohs.loc[two].min()) / m
    intrinsic = pd.DataFrame(intrinsic).T
    drop = {"Co2": ["cell_045", "cell_046", "cell_069", "cell_070"], "Co5": ["cell_087", "cell_088"],
            "Co10": ["cell_071", "cell_072"]}
    groups = {}
    for crate in ("Co10", "Co5", "Co2"):
        dyn = meta[(meta["avg_Crate"] == crate) & (meta["protocol_type"] != "CC")]["cell_name"].tolist()
        nost = meta[(meta["avg_Crate"] == crate) & meta["protocol_variant"].str.contains("no storage")]["cell_name"].tolist()
        bucket = [c for c in dyn if c not in drop[crate]] + nost
        pois = [p for p in meta[meta["cell_name"].isin(bucket)]["protocol_name"].unique() if p not in SINGLETON_PROTOCOLS]
        sub = intrinsic.loc[pois]
        table = all_sohs.loc[bucket]
        cc = meta[(meta["avg_Crate"] == crate) & (meta["protocol_type"] == "CC")
                  & ~meta["protocol_variant"].str.contains("% storage")]["cell_name"].values
        cc_mean = all_sohs.loc[cc].mean()
        bars = np.maximum(table.max() - cc_mean, (table.min() - cc_mean).abs())
        cc_mean.iloc[0] = 1
        bars = bars / cc_mean
        groups[CRATE_LABEL[crate]] = {
            "n_protocols": len(pois), "n_cells": len(bucket),
            **{f"soh{s:.3f}": {"cell_to_cell_mean": float(sub[s].mean()), "cell_to_cell_max": float(sub[s].max()),
                               "max_difference_with_constant_current": float(bars[s])} for s in (0.9, 0.875, 0.85)}}
    means = [g[f"soh{s:.3f}"]["cell_to_cell_mean"] for g in groups.values() for s in (0.9, 0.875, 0.85)]
    halves = [g[f"soh{s:.3f}"]["cell_to_cell_max"] / g[f"soh{s:.3f}"]["max_difference_with_constant_current"]
              for g in groups.values() for s in (0.9, 0.875, 0.85)]
    # Our primary pair differences against the notebook's, pair by pair, at the three Fig. 2c SOHs.
    diffs = []
    for s in (0.9, 0.875, 0.85):
        ours = pt[np.isclose(pt["threshold"], s)].set_index("protocol")["delta_rel"]
        theirs = intrinsic[s].dropna()
        common = ours.index.intersection(theirs.index)
        dd = (ours.loc[common] - theirs.loc[common]).abs()
        diffs.append({"soh": s, "pairs_compared": int(len(common)),
                      "max_abs_difference": float(dd.max()),
                      "pairs_differing_by_more_than_0.001": sorted(dd[dd > 1e-3].index.tolist())})
    return {"groups": groups,
            "paper_statement_average_below_5pct": bool(max(means) < 0.05),
            "largest_group_average": float(max(means)),
            "paper_statement_range_at_most_half_of_protocol_spread": bool(max(halves) <= 0.5),
            "largest_range_over_protocol_spread": float(max(halves)),
            "primary_versus_notebook_pair_differences": diffs}


def eol_crosscheck(data_dir, est, pairs):
    """The paper's own lifetime at 90 % SOH (data/eol_metrics_soh0.9.pkl, 'EFCs (with Diagnostic)')
    against ours, and the pair statistics its values give."""
    path = os.path.join(data_dir, "eol_metrics_soh0.9.pkl")
    if not os.path.exists(path):
        return None
    em = read_pickle(path)
    col = "EFCs (with Diagnostic)"
    common = est.index.intersection(em.index)
    theirs = em.loc[common, col].astype(float)
    theirs = theirs[np.isfinite(theirs) & (theirs > 0)]
    ours = est.loc[theirs.index, REFERENCE].astype(float)
    rel = ((ours - theirs) / theirs).dropna()
    d = []
    for p, a, b in pairs:
        if a in theirs.index and b in theirs.index:
            la, lb = float(theirs.loc[a]), float(theirs.loc[b])
            d.append(abs(la - lb) / ((la + lb) / 2.0))
    return {"cells_compared": int(len(rel)), "median_abs_relative_difference": float(rel.abs().median()),
            "max_abs_relative_difference": float(rel.abs().max()),
            "cells_within_1pct": int((rel.abs() <= 0.01).sum()),
            "pair_statistics_from_the_paper_file": summary(d) if len(d) else None}


def sensitivity(est, pairs, rng):
    """Pair statistics, sigma and resampled worst-of-8 and worst-of-96 at 90 and 85 % SOH, for a
    lifetime table (est needs life_soh0.900, life_soh0.850, protocol_type and avg_Crate)."""
    pt = pair_table(est, pairs)
    out = {"n_pairs": len(pairs)}
    for s in (0.90, 0.85):
        d = deltas(pt, s)
        pool = deviation_pool(d)
        res = {**summary(d), **sigma_estimates(d)}
        res["nonparametric_median_worst_of_8"] = float(np.median(worst_of_n_shortfall(rng.choice(pool, (N_MC, 8)))))
        res["nonparametric_median_worst_of_96"] = float(np.median(worst_of_n_shortfall(rng.choice(pool, (N_MC, 96)))))
        out[f"soh{s:.3f}"] = res
    return out


# ---------------------------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------------------------
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def clean(x):
    """JSON-safe copy: numpy scalars to Python, NaN to None."""
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, (np.bool_, bool)):
        return bool(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        return float(x) if math.isfinite(float(x)) else None
    return x


def make_figure(pt, a3, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

    ink, ink2, muted, grid, axis = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
    c90, c85 = "#2a78d6", "#eb6834"
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": axis,
                         "axes.linewidth": 0.8, "xtick.color": ink2, "ytick.color": ink2,
                         "axes.labelcolor": ink, "text.color": ink, "legend.frameon": False})
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.9), dpi=200)
    for a in ax:
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
        a.grid(True, axis="y", color=grid, linewidth=0.6)
        a.set_axisbelow(True)

    # a: cumulative distribution of the pair difference
    xmax = 0.0
    for s, col in ((0.90, c90), (0.85, c85)):
        d = np.sort(deltas(pt, s)) * 100.0
        if not len(d):
            continue
        xmax = max(xmax, d.max())
        y = np.arange(1, len(d) + 1) / len(d)
        ax[0].step(np.r_[0.0, d], np.r_[0.0, y], where="post", color=col, linewidth=1.6,
                   label=f"SOH {s * 100:g} %: {len(d)} pairs, mean {d.mean():.1f} %, 95th pct {np.percentile(d, 95):.1f} %")
    ax[0].axvline(5.0, color=muted, linewidth=0.8, zorder=1)
    ax[0].text(5.0 + 0.012 * max(xmax, 6.0), 0.33, "5 %: the paper's bound\non the average difference",
               color=ink2, fontsize=7.5, va="bottom")
    ax[0].axhline(0.95, color=muted, linewidth=0.6, zorder=1)
    ax[0].text(max(xmax, 6.0) * 1.06, 0.94, "0.95", color=ink2, fontsize=7.5, va="top", ha="right")
    ax[0].set_xlim(0, max(xmax, 6.0) * 1.08)
    ax[0].set_ylim(0, 1.02)
    ax[0].set_xlabel("Relative lifetime difference between duplicate cells (%)")
    ax[0].set_ylabel("Cumulative fraction of pairs")
    ax[0].set_title("a   Duplicate pairs, EFC to the SOH threshold", loc="left", fontsize=9.5)
    ax[0].legend(fontsize=7.5, loc="lower right", bbox_to_anchor=(1.0, 0.02))

    # b: worst of N by selection alone, at 90 % SOH
    key = "soh0.900"
    if key in a3:
        w = a3[key]["worst_of_N_shortfall"]
        ns = [int(n) for n in w]
        get = lambda k: [100.0 * w[str(n)][k] for n in ns]
        ax[1].plot(ns, get("parametric_p95"), color=muted, linewidth=1.0, linestyle=(0, (4, 3)), label="95th percentile, normal model")
        ax[1].plot(ns, get("parametric_median"), color=muted, linewidth=1.0, linestyle=(0, (1, 2)), label="Median, normal model")
        ax[1].plot(ns, get("nonparametric_p95"), color=c90, linewidth=1.6, marker="s", markersize=5,
                   markerfacecolor="white", markeredgewidth=1.4, label="95th percentile, resampled pairs")
        ax[1].plot(ns, get("nonparametric_median"), color=c90, linewidth=1.6, marker="o", markersize=5,
                   markeredgecolor="white", markeredgewidth=1.0, label="Median, resampled pairs")
        for n, ha, dx in ((8, "left", 6), (96, "right", -7)):
            v = 100.0 * w[str(n)]["nonparametric_median"]
            ax[1].annotate(f"{v:.1f} %", (n, v), textcoords="offset points", xytext=(dx, -12), ha=ha,
                           fontsize=7.5, color=ink)
        ax[1].set_xscale("log")
        ax[1].xaxis.set_major_locator(FixedLocator(ns))
        ax[1].xaxis.set_minor_locator(NullLocator())
        ax[1].xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(round(v))}"))
        ax[1].set_xlim(ns[0] / 1.25, ns[-1] * 1.25)
        top = max(100.0 * max(w[str(n)]["nonparametric_p95"], w[str(n)]["parametric_p95"]) for n in ns)
        ax[1].set_ylim(0, top * 1.08)
        ax[1].set_xlabel("Nominally identical cells in the module, N")
        ax[1].set_ylabel("Worst cell below the module mean (%)")
        ax[1].set_title("b   Worst of N by selection alone, SOH 90 %", loc="left", fontsize=9.5)
        ax[1].legend(fontsize=7.5, loc="lower right", bbox_to_anchor=(1.0, 0.02))
    fig.tight_layout()
    fig.savefig(out_path, facecolor="white")
    plt.close(fig)


# ---------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data", required=True, help="the data folder of the Geslin et al. repository")
    ap.add_argument("--out", default="results")
    ap.add_argument("--upstream-commit", default=None,
                    help="commit of the Geslin et al. repository the data folder comes from (recorded only)")
    ap.add_argument("--allow-unverified-inputs", action="store_true",
                    help="run on files whose SHA-256 differs from the released ones (used by tests/synthetic_check.py)")
    args = ap.parse_args()
    hashes = verify_inputs(args.data, args.allow_unverified_inputs)
    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(SEED)

    meta, raw, traj = load(args.data, monotone_filter)
    groups, pairs, singletons = make_pairs(meta)
    est = build_estimators(meta, traj)
    est.to_csv(os.path.join(args.out, "per_cell_estimators.csv"))
    pt = pair_table(est, pairs)
    pt.to_csv(os.path.join(args.out, "pair_differences.csv"), index=False)

    results = {
        "n_cells": int(len(meta)), "n_protocols": int(len(groups)), "n_pairs": int(len(pairs)),
        "singletons": [{"cell": c, "protocol": p, "note": "duplicate failed; excluded from pair statistics by construction"}
                       for c, p in singletons.items()],
        "processing": {"filter": "monotone (every point that a later point exceeds is dropped)",
                       "lifetime": "EFC at SOH x first retained C/2 capacity, linear interpolation; NaN if not reached",
                       "thresholds": THRESHOLDS, "seed": SEED},
    }
    results["A1_pair_relative_difference"] = a1_summary(pt)
    a2 = a2_sigma(pt, rng)
    results["A2_cell_noise_sigma"] = a2
    results["A2_tail_diagnostics"] = {f"soh{s:.3f}": tail_diagnostics(deltas(pt, s), rng) for s in (0.90, 0.85)}
    a3 = a3_worst(pt, a2, rng)
    results["A3_selection_only_worst_of_N"] = a3
    results["headline_percentile_bootstrap_95"] = headline_bootstrap(deltas(pt, 0.90), deltas(pt, 0.85), rng)
    a4, agree, rr = a4_estimators(est, pairs)
    agree.to_csv(os.path.join(args.out, "worse_cell_agreement_matrix.csv"))
    rr.to_csv(os.path.join(args.out, "pair_ratio_ranges.csv"), index=False)
    results["A4_estimator_sensitivity"] = a4
    results["A5_variance_split"] = a5_variance(est)
    nb = notebook_lifetimes(meta, raw)
    results["R_reconciliation"] = {"paper_fig2c": paper_fig2c(meta, nb, pt),
                                   "paper_eol_metrics_soh0.9": eol_crosscheck(args.data, est, pairs)}
    est_nb = pd.DataFrame({f"life_soh{s:.3f}": nb[float(np.round(s, 3))] for s in THRESHOLDS})
    est_nb = est_nb.join(meta.set_index("cell_name")[["protocol_name", "protocol_type", "avg_Crate"]])
    results["S1_paper_notebook_processing"] = sensitivity(est_nb, pairs, rng)
    _, pairs_x, _ = make_pairs(meta, exclude=C16_DRIVE_PROTOCOLS)
    results["S2_without_c16_drive_protocols"] = sensitivity(est, pairs_x, rng)
    _, pairs_y, _ = make_pairs(meta, exclude=AUTHORS_ANALYSIS_SET_EXCLUDED)
    results["S3_authors_analysis_set_without_c16_drive_and_periodic_c_c2"] = sensitivity(est, pairs_y, rng)
    results["sanity"] = {
        "n_pairs_at_soh90": int(len(deltas(pt, 0.90))),
        "n_pairs_at_soh85": int(len(deltas(pt, 0.85))),
        "mean_pair_difference_at_soh90": float(deltas(pt, 0.90).mean()),
        "within_pair_variance_share_at_soh90": results["A5_variance_split"].get("soh0.900", {}).get("within_pair_variance_share"),
        "paper_fig2c_average_below_5pct_reproduced": results["R_reconciliation"]["paper_fig2c"]["paper_statement_average_below_5pct"],
    }
    try:
        import matplotlib
        mpl_version = matplotlib.__version__
    except Exception:
        mpl_version = None
    results["provenance"] = {
        "inputs_sha256": hashes,
        "inputs_match_the_released_files": all(hashes.get(f) == h for f, h in EXPECTED_SHA256.items()),
        "upstream_commit": args.upstream_commit,
        "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
        "matplotlib": mpl_version,
    }
    with open(os.path.join(args.out, "results.json"), "w") as f:
        json.dump(clean(results), f, indent=2, allow_nan=False)

    try:
        make_figure(pt, a3, os.path.join(args.out, "lifetime_floor.png"))
    except Exception as ex:  # the figure is optional; the numbers are not
        with open(os.path.join(args.out, "figure_error.txt"), "w") as f:
            f.write(repr(ex))

    a1 = results["A1_pair_relative_difference"]
    w90 = a3.get("soh0.900", {}).get("worst_of_N_shortfall", {})
    print(f"pairs: {len(pairs)} ({len(deltas(pt, 0.90))} at SOH 90 %, {len(deltas(pt, 0.85))} at SOH 85 %); "
          f"singletons: {', '.join(singletons)}")
    for k in ("soh0.900", "soh0.850"):
        if k in a1:
            print(f"{k}: mean {100 * a1[k]['mean']:.2f} %, median {100 * a1[k]['median']:.2f} %, "
                  f"95th pct {100 * a1[k]['p95']:.2f} %, max {100 * a1[k]['max']:.2f} %")
    if w90:
        print(f"worst of 8 / 96 at SOH 90 % (resampled median): {100 * w90['8']['nonparametric_median']:.2f} % / "
              f"{100 * w90['96']['nonparametric_median']:.2f} %")
    print("done; results in", args.out)


if __name__ == "__main__":
    main()
