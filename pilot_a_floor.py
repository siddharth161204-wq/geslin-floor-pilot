#!/usr/bin/env python3
"""
Pilot A: the cell-to-cell floor in the Geslin et al. (Nature Energy 2025) dataset.

What it computes, from the processed per-cell diagnostics released with the paper
(github.com/geslina/dynamic_cycling_Nature_Energy_2024, data/):

  A1  Duplicate-pair relative lifetime difference at several SOH thresholds:
      delta = |EFC_a - EFC_b| / mean(EFC_a, EFC_b), with the full distribution
      (mean, median, 75th, 90th, 95th percentiles, max), overall and by protocol type.
  A2  The per-cell lifetime noise sigma implied by the pair differences, with a
      bootstrap confidence interval (pairs resampled).
  A3  The selection-only "worst of N" shortfall: in a module of N nominally identical
      cells under one protocol at one temperature, how far below the mean the worst
      cell's lifetime is expected to fall by selection alone (median and 95th
      percentile), for N = 2 ... 96, parametric (normal) and nonparametric.
  A4  Estimator sensitivity: which duplicate is the "worse" cell, and the pair ratio,
      under several lifetime / rate estimators; fraction of pairs whose worse cell
      flips between estimators; per-pair range of the ratio across estimators.
  A5  Variance split: share of the between-cell lifetime variance at each threshold
      that lies between protocols versus within duplicate pairs.

Inputs are read only; nothing in the Geslin repository is modified.
Processing of the trajectories follows the paper's own notebook (figure2.ipynb):
the capacity curve is made strictly decreasing by dropping invalid points, then
EFC at each SOH is found by linear interpolation; cells that never reach a
threshold are NaN at that threshold.

Usage:
  python3 pilot_a_floor.py --data /path/to/dynamic_cycling_Nature_Energy_2024/data --out results
"""
import argparse, json, os
import numpy as np
import pandas as pd

RNG = np.random.default_rng(20261007)

CAPA = "C/2 discharge CC Capacity [normalized]"
EFC = "EFCs"
THRESHOLDS = [0.95, 0.925, 0.90, 0.875, 0.85]
N_LIST = [2, 4, 8, 12, 24, 48, 96]
N_BOOT = 5000
N_MC = 20000


def remove_invalid_points(cycle, capacity):
    """Paper's own filter: drop points that break strict monotonic decrease."""
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


def efc_at_soh(efc, capa, thresholds):
    """Linear interpolation of EFC at SOH*capa[0]; NaN if not reached."""
    capa0 = capa[0]
    out = []
    for s in thresholds:
        target = s * capa0
        if capa.min() > target:
            out.append(np.nan)
            continue
        out.append(float(np.interp(target, np.flip(capa), np.flip(efc), left=np.nan, right=np.nan)))
    return np.array(out)


def load(data_dir):
    meta = pd.read_pickle(os.path.join(data_dir, "metadata.pkl"))
    dfa = pd.read_pickle(os.path.join(data_dir, "diagnostic_features_all.pkl"))
    cells = list(meta["cell_name"])
    traj = {}
    for cell in cells:
        sub = dfa.loc[cell]
        e = np.array(sub[EFC], dtype=float)
        q = np.array(sub[CAPA].values, dtype=float)
        e, q = remove_invalid_points(e, q)
        traj[cell] = (e, q)
    return meta, traj


def estimators(e, q):
    """Several ways to say how fast a cell degrades. Lower 'life' = worse; higher 'rate' = worse."""
    out = {}
    life = efc_at_soh(e, q, THRESHOLDS)
    for s, v in zip(THRESHOLDS, life):
        out[f"life_soh{s:.3f}"] = v                      # lifetime estimators (higher = better)
    q0 = q[0]
    soh = q / q0
    # average fade rate up to the last diagnostic (per 100 EFC; higher = worse)
    out["rate_avg_last"] = -100.0 * (soh[-1] - 1.0) / e[-1] if e[-1] > 0 else np.nan
    # linear slope over the first k diagnostics (higher = worse)
    for k in (4, 6, 8):
        if len(e) >= k:
            A = np.vstack([e[:k], np.ones(k)]).T
            m, _ = np.linalg.lstsq(A, soh[:k], rcond=None)[0]
            out[f"rate_slope_first{k}"] = -100.0 * m
        else:
            out[f"rate_slope_first{k}"] = np.nan
    # slope over the whole trajectory
    A = np.vstack([e, np.ones(len(e))]).T
    m, _ = np.linalg.lstsq(A, soh, rcond=None)[0]
    out["rate_slope_all"] = -100.0 * m
    # fade at a fixed EFC (interpolated; higher fade = worse)
    for efc_fixed in (200.0, 400.0):
        if e[-1] >= efc_fixed:
            out[f"fade_at_{int(efc_fixed)}efc"] = 100.0 * (1.0 - float(np.interp(efc_fixed, e, soh)))
        else:
            out[f"fade_at_{int(efc_fixed)}efc"] = np.nan
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", default="results")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    meta, traj = load(args.data)
    cells = list(meta["cell_name"])
    est = pd.DataFrame({c: estimators(*traj[c]) for c in cells}).T
    est.index.name = "cell"
    est = est.join(meta.set_index("cell_name")[["protocol_name", "protocol_type", "avg_Crate"]])
    est.to_csv(os.path.join(args.out, "per_cell_estimators.csv"))

    # ---- pairs ----
    groups = meta.groupby("protocol_name")["cell_name"].apply(list)
    pairs = [(p, v[0], v[1]) for p, v in groups.items() if len(v) == 2]
    results = {"n_cells": len(cells), "n_protocols": int(len(groups)), "n_pairs": len(pairs),
               "singletons": [v[0] for v in groups if len(v) == 1]}

    # A1: pair relative difference per threshold
    a1 = {}
    pair_rows = []
    for s in THRESHOLDS:
        col = f"life_soh{s:.3f}"
        deltas, types, crates = [], [], []
        for p, a, b in pairs:
            la, lb = est.loc[a, col], est.loc[b, col]
            if np.isnan(la) or np.isnan(lb):
                continue
            d = abs(la - lb) / ((la + lb) / 2.0)
            deltas.append(d); types.append(est.loc[a, "protocol_type"]); crates.append(est.loc[a, "avg_Crate"])
            pair_rows.append({"threshold": s, "protocol": p, "cell_a": a, "cell_b": b, "life_a": la, "life_b": lb, "delta_rel": d,
                              "protocol_type": types[-1], "avg_Crate": crates[-1]})
        deltas = np.array(deltas)
        if len(deltas) == 0:
            continue
        summ = {"n_pairs": int(len(deltas)), "mean": float(deltas.mean()), "median": float(np.median(deltas)),
                "p75": float(np.percentile(deltas, 75)), "p90": float(np.percentile(deltas, 90)),
                "p95": float(np.percentile(deltas, 95)), "max": float(deltas.max())}
        by_type = {}
        for t in sorted(set(types)):
            dd = deltas[np.array(types) == t]
            by_type[t] = {"n": int(len(dd)), "mean": float(dd.mean()), "median": float(np.median(dd)), "p95": float(np.percentile(dd, 95)) if len(dd) >= 2 else None}
        by_crate = {}
        for t in sorted(set(crates)):
            dd = deltas[np.array(crates) == t]
            by_crate[t] = {"n": int(len(dd)), "mean": float(dd.mean()), "median": float(np.median(dd))}
        summ["by_protocol_type"] = by_type
        summ["by_crate"] = by_crate
        a1[f"soh{s:.3f}"] = summ
    results["A1_pair_relative_difference"] = a1
    pd.DataFrame(pair_rows).to_csv(os.path.join(args.out, "pair_differences.csv"), index=False)

    # A2: per-cell noise sigma from pair differences, with bootstrap CI
    a2 = {}
    for s in THRESHOLDS:
        col = f"life_soh{s:.3f}"
        d = np.array([r["delta_rel"] for r in pair_rows if r["threshold"] == s])
        if len(d) < 5:
            continue
        # E|eps_a - eps_b| = 2 sigma / sqrt(pi) for normal eps; median|.| = 0.6745 * sqrt(2) sigma
        sig_mean = d.mean() * np.sqrt(np.pi) / 2.0
        sig_med = np.median(d) / (0.6745 * np.sqrt(2.0))
        boots = []
        for _ in range(N_BOOT):
            db = RNG.choice(d, size=len(d), replace=True)
            boots.append(db.mean() * np.sqrt(np.pi) / 2.0)
        boots = np.array(boots)
        a2[f"soh{s:.3f}"] = {"sigma_from_mean": float(sig_mean), "sigma_from_median": float(sig_med),
                             "sigma_ci95": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))], "n_pairs": int(len(d))}
    results["A2_cell_noise_sigma"] = a2

    # A3: worst-of-N selection shortfall (parametric and nonparametric), at 90 % and 85 % SOH
    a3 = {}
    for s in (0.90, 0.85):
        key = f"soh{s:.3f}"
        if key not in a2:
            continue
        sig = a2[key]["sigma_from_mean"]
        d = np.array([r["delta_rel"] for r in pair_rows if r["threshold"] == s])
        # nonparametric eps pool: each cell's deviation from its pair mean is +-delta/2; pair mean removes half the variance,
        # so scale by sqrt(2) to recover a single-cell deviation
        eps_pool = np.concatenate([d / 2.0, -d / 2.0]) * np.sqrt(2.0)
        per_n = {}
        for n in N_LIST:
            eps_par = RNG.normal(0.0, sig, size=(N_MC, n))
            short_par = -(eps_par.min(axis=1) - eps_par.mean(axis=1))  # (mean - min), relative
            eps_np = RNG.choice(eps_pool, size=(N_MC, n), replace=True)
            short_np = -(eps_np.min(axis=1) - eps_np.mean(axis=1))
            per_n[str(n)] = {"parametric_median": float(np.median(short_par)), "parametric_p95": float(np.percentile(short_par, 95)),
                             "nonparametric_median": float(np.median(short_np)), "nonparametric_p95": float(np.percentile(short_np, 95))}
        a3[key] = {"sigma_used": float(sig), "worst_of_N_shortfall": per_n}
    results["A3_selection_only_worst_of_N"] = a3

    # A4: estimator sensitivity on pairs
    life_cols = [f"life_soh{s:.3f}" for s in THRESHOLDS]
    rate_cols = [c for c in est.columns if c.startswith("rate_") or c.startswith("fade_")]
    worse = {}
    ratio_range = []
    for p, a, b in pairs:
        w = {}
        ratios = {}
        for c in life_cols:
            la, lb = est.loc[a, c], est.loc[b, c]
            if np.isnan(la) or np.isnan(lb):
                continue
            w[c] = a if la < lb else b
            ratios[c] = max(la, lb) / min(la, lb)
        for c in rate_cols:
            ra, rb = est.loc[a, c], est.loc[b, c]
            if np.isnan(ra) or np.isnan(rb) or min(ra, rb) <= 0:
                continue
            w[c] = a if ra > rb else b
            ratios[c] = max(ra, rb) / min(ra, rb)
        worse[p] = w
        if len(ratios) >= 3:
            ratio_range.append({"protocol": p, "n_estimators": len(ratios), "ratio_min": float(min(ratios.values())),
                                "ratio_max": float(max(ratios.values()))})
    est_names = life_cols + rate_cols
    agree = pd.DataFrame(index=est_names, columns=est_names, dtype=float)
    for c1 in est_names:
        for c2 in est_names:
            both = [p for p in worse if c1 in worse[p] and c2 in worse[p]]
            agree.loc[c1, c2] = np.mean([worse[p][c1] == worse[p][c2] for p in both]) if both else np.nan
    agree.to_csv(os.path.join(args.out, "worse_cell_agreement_matrix.csv"))
    flips = []
    for p, w in worse.items():
        ids = set(w.values())
        if len(w) >= 3:
            flips.append(len(ids) > 1)
    rr = pd.DataFrame(ratio_range)
    results["A4_estimator_sensitivity"] = {
        "n_pairs_with_3plus_estimators": int(len(flips)),
        "fraction_of_pairs_whose_worse_cell_depends_on_estimator": float(np.mean(flips)) if flips else None,
        "pair_ratio_range_median_min": float(rr["ratio_min"].median()) if len(rr) else None,
        "pair_ratio_range_median_max": float(rr["ratio_max"].median()) if len(rr) else None,
        "pair_ratio_range_overall": [float(rr["ratio_min"].min()), float(rr["ratio_max"].max())] if len(rr) else None,
        "agreement_life90_vs_life85": float(agree.loc["life_soh0.900", "life_soh0.850"]),
        "agreement_life90_vs_slope_first6": float(agree.loc["life_soh0.900", "rate_slope_first6"]),
    }
    rr.to_csv(os.path.join(args.out, "pair_ratio_ranges.csv"), index=False)

    # A5: variance split between protocols and within pairs
    a5 = {}
    for s in THRESHOLDS:
        col = f"life_soh{s:.3f}"
        vals = est[[col, "protocol_name"]].dropna()
        g = vals.groupby("protocol_name")[col]
        twos = [v.values for _, v in g if len(v) == 2]
        if len(twos) < 5:
            continue
        within = np.mean([np.var(v, ddof=1) for v in twos])
        allv = np.concatenate(twos)
        total = np.var(allv, ddof=1)
        a5[f"soh{s:.3f}"] = {"within_pair_variance_share": float(within / total), "between_protocol_variance_share": float(1 - within / total),
                             "n_pairs": len(twos)}
    results["A5_variance_split"] = a5

    with open(os.path.join(args.out, "results.json"), "w") as f:
        json.dump(results, f, indent=2)

    # ---- figure ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(9, 3.6), dpi=150)
        for s, ls in ((0.90, "-"), (0.85, "--")):
            d = np.sort([r["delta_rel"] for r in pair_rows if r["threshold"] == s])
            if len(d):
                ax[0].step(100 * d, np.arange(1, len(d) + 1) / len(d), ls, label=f"SOH {int(s*100)} % (n = {len(d)} pairs)")
        ax[0].axvline(5, color="gray", lw=0.8); ax[0].text(5.2, 0.05, "5 %", color="gray", fontsize=8)
        ax[0].set_xlabel("Duplicate relative lifetime difference (%)"); ax[0].set_ylabel("Cumulative fraction of pairs"); ax[0].legend(fontsize=8)
        if "soh0.900" in a3:
            ns = [int(n) for n in a3["soh0.900"]["worst_of_N_shortfall"]]
            med = [100 * a3["soh0.900"]["worst_of_N_shortfall"][str(n)]["nonparametric_median"] for n in ns]
            p95 = [100 * a3["soh0.900"]["worst_of_N_shortfall"][str(n)]["nonparametric_p95"] for n in ns]
            ax[1].plot(ns, med, "o-", label="median (nonparametric)")
            ax[1].plot(ns, p95, "s--", label="95th percentile")
            ax[1].set_xscale("log"); ax[1].set_xlabel("Cells in module, N"); ax[1].set_ylabel("Worst cell below mean by selection alone (%)")
            ax[1].legend(fontsize=8); ax[1].set_title("Floor at SOH 90 %", fontsize=9)
        fig.tight_layout(); fig.savefig(os.path.join(args.out, "pilot_a_floor.png")); plt.close(fig)
    except Exception as ex:  # figure is optional
        with open(os.path.join(args.out, "figure_error.txt"), "w") as f:
            f.write(str(ex))

    print("done; results in", args.out)


if __name__ == "__main__":
    main()