"""Independent recomputation of pair-variability numbers from the Geslin et al.
(Nature Energy 2025) dynamic-cycling dataset. Written from the task definitions only.

Run:  python3 -I verification/recompute.py <path to dynamic_cycling_Nature_Energy_2024/data>
Writes verification/recompute.json next to this file.
"""
import hashlib
import json
import math
import sys

import numpy as np
import pandas as pd

import os  # noqa: E402
DATA = (sys.argv[1] if len(sys.argv) > 1 else "dynamic_cycling_Nature_Energy_2024/data").rstrip("/") + "/"
OUT = os.path.dirname(os.path.abspath(__file__)) + "/"
EXPECTED = {
    "metadata.pkl": "7579f48cf3767552154fba3ad9acb048c3166088c2a46fadb06b2a221b3a1e1f",
    "diagnostic_features_all.pkl": "906e54db95cf9a3f294b8fa0ba29d43b5850a7f1a228cdd9e37cd2bcdfe1bc79",
}
QCOL = "C/2 discharge CC Capacity [normalized]"
ECOL = "EFCs"
LIFE_S = [0.95, 0.925, 0.90, 0.875, 0.85]
NAN = float("nan")


# ----------------------------------------------------------------------------- load
def verify_and_load():
    for fn, h in EXPECTED.items():
        with open(DATA + fn, "rb") as f:
            got = hashlib.sha256(f.read()).hexdigest()
        if got != h:
            sys.exit(f"SHA-256 mismatch for {fn}: {got}")
    meta = pd.read_pickle(DATA + "metadata.pkl")
    diag = pd.read_pickle(DATA + "diagnostic_features_all.pkl")
    return meta, diag


def cell_series(diag):
    """Per-cell (q, e) arrays in stored row order."""
    lvl0 = diag.index.get_level_values(0)
    qv = diag[QCOL].to_numpy(dtype=float)
    ev = diag[ECOL].to_numpy(dtype=float)
    out = {}
    for cell in pd.unique(lvl0):
        m = np.asarray(lvl0 == cell)
        out[cell] = (qv[m].copy(), ev[m].copy())
    return out


# ----------------------------------------------------------------------------- filters
def filter_primary(q, e):
    """Keep point i only if q[i] >= q[j] for every j > i (drop any point exceeded by ANY later point)."""
    n = len(q)
    keep = np.zeros(n, dtype=bool)
    later_max = -np.inf
    for i in range(n - 1, -1, -1):
        if q[i] >= later_max:
            keep[i] = True
        later_max = max(later_max, q[i])
    # brute-force confirmation of the definition
    brute = np.array([bool(np.all(q[i] >= q[i + 1:])) for i in range(n)])
    assert np.array_equal(keep, brute)
    return q[keep], e[keep]


def filter_runmin_start(q, e):
    """Sensitivity variant only: keep point i only if q[i] <= q[j] for every j < i."""
    n = len(q)
    keep = np.array([bool(np.all(q[i] <= q[:i])) for i in range(n)])
    return q[keep], e[keep]


def filter_notebook(q, e):
    """figure2.ipynb remove_invalid_points: original point i-1 is dropped whenever q[i] > q[i-1]
    (adjacent comparison on the ORIGINAL series; the backward loop never shifts lower indices)."""
    keep = np.append(~(q[1:] > q[:-1]), True)
    return q[keep], e[keep]


# ----------------------------------------------------------------------------- lifetime
def lifetime(qf, ef, s):
    """EFC at which the filtered curve first reaches s*qf[0]; NaN if min(qf) > target."""
    target = s * qf[0]
    if np.min(qf) > target:
        return NAN
    k = int(np.argmax(qf <= target))  # first point at or below target
    if k == 0:
        return float(ef[0])
    q0, q1 = qf[k - 1], qf[k]
    e0, e1 = ef[k - 1], ef[k]
    return float(e0 + (target - q0) * (e1 - e0) / (q1 - q0))


def pair_diff(la, lb):
    if math.isnan(la) or math.isnan(lb):
        return NAN
    return abs(la - lb) / ((la + lb) / 2.0)


# ----------------------------------------------------------------------------- estimators
def ls_slope(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    xm = x.mean()
    ym = y.mean()
    return float(np.sum((x - xm) * (y - ym)) / np.sum((x - xm) ** 2))


ESTIMATORS = [
    # name, direction ('low' = lower is worse, 'high' = higher is worse), needs_positive
    ("L0.95", "low", False), ("L0.925", "low", False), ("L0.90", "low", False),
    ("L0.875", "low", False), ("L0.85", "low", False),
    ("avg_fade_rate", "high", True),
    ("neg_slope_first4", "high", True), ("neg_slope_first6", "high", True),
    ("neg_slope_first8", "high", True), ("neg_slope_all", "high", True),
    ("fade_at_200", "high", True), ("fade_at_400", "high", True),
]


def estimators_for_cell(qf, ef):
    soh = qf / qf[0]
    v = {}
    for s, name in zip(LIFE_S, ["L0.95", "L0.925", "L0.90", "L0.875", "L0.85"]):
        v[name] = lifetime(qf, ef, s)
    v["avg_fade_rate"] = float((1.0 - soh[-1]) / (ef[-1] - ef[0])) if ef[-1] != ef[0] else NAN
    for k in (4, 6, 8):
        v[f"neg_slope_first{k}"] = -ls_slope(ef[:k], soh[:k]) if len(ef) >= k else NAN
    v["neg_slope_all"] = -ls_slope(ef, soh) if len(ef) >= 2 else NAN
    for x in (200.0, 400.0):
        if ef[0] <= x <= ef[-1]:
            v[f"fade_at_{int(x)}"] = float(1.0 - np.interp(x, ef, soh))
        else:
            v[f"fade_at_{int(x)}"] = NAN
    return v


def pair_verdicts(ca, cb, va, vb):
    """Return (n_defined_both, {estimator: worse_cell}) following the task rules."""
    verdicts = {}
    n_def = 0
    for name, direction, needs_pos in ESTIMATORS:
        a, b = va[name], vb[name]
        if math.isnan(a) or math.isnan(b):
            continue
        n_def += 1
        if a == b:
            continue
        if needs_pos and not (a > 0 and b > 0):
            continue
        if direction == "low":
            verdicts[name] = ca if a < b else cb
        else:
            verdicts[name] = ca if a > b else cb
    return n_def, verdicts


# ----------------------------------------------------------------------------- worst-of-N
def worst_of_n_median(d, N, n_modules, seed, chunk=50_000, return_samples=False):
    pool = np.concatenate([d / 2.0, -d / 2.0]) * np.sqrt(2.0)
    rng = np.random.default_rng(seed)
    shortfall = np.empty(n_modules)
    for start in range(0, n_modules, chunk):
        m = min(chunk, n_modules - start)
        idx = rng.integers(0, pool.size, size=(m, N))
        life = 1.0 + pool[idx]
        mean = life.mean(axis=1)
        shortfall[start:start + m] = (mean - life.min(axis=1)) / mean
    if return_samples:
        return float(np.median(shortfall)), pool, shortfall
    return float(np.median(shortfall)), pool


# ----------------------------------------------------------------------------- figure-2c reimplementation
def fig2c_reimplementation(meta, series):
    nb_cells = ["cell_" + "{:03d}".format(c) for c in range(3, 97) if c not in (18, 83)]
    assert sorted(nb_cells) == sorted(meta["cell_name"].tolist())
    SOHs = np.linspace(1, 0.85, 31)
    cols = np.round(SOHs, 3)
    eol = {}
    nonmono = []
    for cell in nb_cells:
        q, e = series[cell]
        qf, ef = filter_notebook(q, e)
        if np.any(np.diff(qf) > 0):
            nonmono.append(cell)
        # same call, same argument order (x descending) as the notebook
        v = np.round(np.interp(SOHs * qf[0], np.flip(qf), np.flip(ef), left=np.nan, right=np.nan), 1)
        v[0] = 0
        eol[cell] = v
    # protocol (max-min)/mean with pandas skipna semantics written out explicitly
    single = {"CC_(no_storage)_Co2", "Synthetic_2c_Co5"}
    prots = [p for p in pd.unique(meta["protocol_name"]) if p not in single]
    intrinsic = {}
    one_nan = {}
    for p in prots:
        two = meta.loc[meta["protocol_name"] == p, "cell_name"].tolist()
        assert len(two) == 2
        a, b = eol[two[0]], eol[two[1]]
        row = np.empty(len(SOHs))
        for j in range(len(SOHs)):
            vals = [x for x in (a[j], b[j]) if not np.isnan(x)]
            if len(vals) == 0:
                row[j] = np.nan
            else:
                mean = 1.0 if j == 0 else sum(vals) / len(vals)
                row[j] = (max(vals) - min(vals)) / mean
            if len(vals) == 1:
                one_nan.setdefault(p, []).append(float(cols[j]))
        intrinsic[p] = row
    # buckets as in the notebook
    crates = ["Co10", "Co5", "Co2"]
    excl = {"Co2": ["cell_045", "cell_046", "cell_069", "cell_070"],
            "Co5": ["cell_087", "cell_088"], "Co10": ["cell_071", "cell_072"]}
    groups = {}
    for cr in crates:
        dyn = meta.loc[(meta["avg_Crate"] == cr) & (meta["protocol_type"] != "CC"), "cell_name"].tolist()
        dyn = [c for c in dyn if c not in excl[cr]]
        nost = meta.loc[(meta["avg_Crate"] == cr) & (meta["protocol_variant"].str.contains("no storage")), "cell_name"].tolist()
        b3 = dyn + nost
        pr = pd.unique(meta.loc[meta["cell_name"].isin(b3), "protocol_name"])
        pr = [p for p in pr if p not in single]
        res = {"n_protocols": len(pr), "protocols": pr, "bucket_cells": b3}
        for soh in (0.90, 0.875, 0.85):
            j = int(np.where(cols == soh)[0][0])
            vals = np.array([intrinsic[p][j] for p in pr])
            ok = ~np.isnan(vals)
            res[f"{soh}"] = {
                "mean": float(vals[ok].mean()),
                "n_nonnan": int(ok.sum()),
                "n_protocols_with_one_cell_nan_counted_as_0": int(sum(1 for p in pr if soh in one_nan.get(p, []))),
                "max": float(vals[ok].max()),
            }
        groups[f"C/{cr[2:]}"] = res
    return groups, nonmono, one_nan


# ----------------------------------------------------------------------------- main analysis per filter
def analyse(meta, series, filt, with_estimators=True):
    filtered = {c: filt(*series[c]) for c in meta["cell_name"]}
    groups = meta.groupby("protocol_name", sort=False)["cell_name"].apply(list)
    pairs = [(p, cs) for p, cs in groups.items() if len(cs) == 2]
    singles = [p for p, cs in groups.items() if len(cs) != 2]
    life = {c: {s: lifetime(*filtered[c], s) for s in LIFE_S} for c in meta["cell_name"]}
    res = {"n_pairs": len(pairs), "single_cell_protocols": singles, "per_threshold": {}}
    for s in (0.90, 0.85, 0.95, 0.925, 0.875):
        rows = []
        for p, (ca, cb) in pairs:
            d = pair_diff(life[ca][s], life[cb][s])
            rows.append((p, ca, cb, life[ca][s], life[cb][s], d))
        valid = [r for r in rows if not math.isnan(r[5])]
        d = np.array([r[5] for r in valid])
        imax = int(np.argmax(d))
        top5 = sorted(valid, key=lambda r: -r[5])[:5]
        la = np.array([r[3] for r in valid])
        lb = np.array([r[4] for r in valid])
        W = float(np.mean([np.var([x, y], ddof=1) for x, y in zip(la, lb)]))
        T = float(np.var(np.concatenate([la, lb]), ddof=1))
        res["per_threshold"][s] = {
            "n_pairs_both": len(valid),
            "pairs_missing": [(r[0], r[1], r[2], r[3], r[4]) for r in rows if math.isnan(r[5])],
            "d": d,
            "mean": float(d.mean()), "median": float(np.median(d)),
            "p95": float(np.percentile(d, 95)), "max": float(d[imax]),
            "max_protocol": valid[imax][0],
            "sigma": float(d.mean() * math.sqrt(math.pi) / 2.0),
            "within_pair_share": W / T, "W": W, "T": T,
            "top5": [(r[0], r[1], r[2], r[3], r[4], r[5]) for r in top5],
            "rows": rows,
        }
    if with_estimators:
        est = {c: estimators_for_cell(*filtered[c]) for c in meta["cell_name"]}
        dep_rows = []
        for p, (ca, cb) in pairs:
            n_def, verd = pair_verdicts(ca, cb, est[ca], est[cb])
            names = set(verd.values())
            dep_main = (len(verd) >= 3) and (len(names) > 1)       # reading 1: >=3 estimators naming a worse cell
            dep_alt = (n_def >= 3) and (len(names) > 1)            # reading 2: >=3 estimators defined for both cells
            dep_rows.append({"protocol": p, "cells": [ca, cb], "n_defined_both": n_def,
                             "n_verdicts": len(verd), "verdicts": verd,
                             "depends": dep_main, "depends_alt_reading": dep_alt})
        res["estimators"] = {
            "per_pair": dep_rows,
            "n_depends": sum(r["depends"] for r in dep_rows),
            "n_depends_alt_reading": sum(r["depends_alt_reading"] for r in dep_rows),
            "n_pairs": len(dep_rows),
            "per_cell": est,
        }
    return res, filtered


def pct3(x):
    return float(f"{100.0 * x:.3f}")


def main():
    meta, diag = verify_and_load()
    series = cell_series(diag)
    assert set(series) == set(meta["cell_name"])

    prim, filtered = analyse(meta, series, filter_primary, with_estimators=True)
    alt_b, _ = analyse(meta, series, filter_runmin_start, with_estimators=True)
    alt_nb, _ = analyse(meta, series, filter_notebook, with_estimators=True)
    assert prim["n_pairs"] == 45

    # worst-of-N at 0.90
    d90 = prim["per_threshold"][0.90]["d"]
    SEED = 20261007
    NMOD = 10_000_000
    NBATCH = 10
    w8, pool, s8 = worst_of_n_median(d90, 8, NMOD, SEED, chunk=100_000, return_samples=True)
    w96, _, s96 = worst_of_n_median(d90, 96, NMOD, SEED + 1, chunk=100_000, return_samples=True)
    # batch-means standard error of the median (10 batches of 1,000,000 modules)
    se = {}
    for N, s in ((8, s8), (96, s96)):
        bm = np.array([np.median(b) for b in np.split(s, NBATCH)])
        se[N] = float(bm.std(ddof=1) / math.sqrt(NBATCH))
    del s8, s96
    spread = {8: [], 96: []}
    for extra in (1, 2, 3, 4, 5):
        spread[8].append(worst_of_n_median(d90, 8, 1_000_000, 1000 + extra)[0])
        spread[96].append(worst_of_n_median(d90, 96, 1_000_000, 2000 + extra)[0])

    fig2c, nb_nonmono, nb_one_nan = fig2c_reimplementation(meta, series)

    P = prim["per_threshold"]
    report = {
        "environment": {"numpy": np.__version__, "pandas": pd.__version__, "python": sys.version.split()[0]},
        "data_sha256_verified": EXPECTED,
        "n_pairs": prim["n_pairs"],
        "single_cell_protocols": prim["single_cell_protocols"],
        "item1_pairs_with_both_lifetimes": {"0.90": P[0.90]["n_pairs_both"], "0.85": P[0.85]["n_pairs_both"],
                                            "pairs_missing_at_0.85": P[0.85]["pairs_missing"],
                                            "pairs_missing_at_0.90": P[0.90]["pairs_missing"]},
        "item2_pair_difference_percent": {
            str(s): {"mean": pct3(P[s]["mean"]), "median": pct3(P[s]["median"]),
                     "p95_numpy_linear": pct3(P[s]["p95"]), "max": pct3(P[s]["max"]),
                     "max_protocol": P[s]["max_protocol"],
                     "unrounded_fraction": {"mean": P[s]["mean"], "median": P[s]["median"],
                                            "p95": P[s]["p95"], "max": P[s]["max"]}}
            for s in (0.90, 0.85)},
        "item3_sigma": {str(s): {"fraction": P[s]["sigma"], "percent_3dp": pct3(P[s]["sigma"])} for s in (0.90, 0.85)},
        "item4_worst_of_N_median_shortfall_at_0.90": {
            "seed_N8": SEED, "seed_N96": SEED + 1, "n_modules_each": NMOD, "rng": "numpy.random.default_rng (PCG64), integers() index draws",
            "pool_size": int(pool.size),
            "N8": {"fraction": w8, "percent_3dp": pct3(w8), "mc_standard_error_percent": 100.0 * se[8]},
            "N96": {"fraction": w96, "percent_3dp": pct3(w96), "mc_standard_error_percent": 100.0 * se[96]},
            "mc_check_5_other_seeds_1M_modules_each_percent": {"seeds_N8": [1001, 1002, 1003, 1004, 1005],
                                                               "seeds_N96": [2001, 2002, 2003, 2004, 2005],
                                                               "N8": [pct3(x) for x in spread[8]],
                                                               "N96": [pct3(x) for x in spread[96]]},
        },
        "item5_within_pair_variance_share_at_0.90": {"share": P[0.90]["within_pair_share"], "W": P[0.90]["W"],
                                                      "T": P[0.90]["T"], "n_pairs": P[0.90]["n_pairs_both"]},
        "item6_estimator_dependence": {
            "n_depends": prim["estimators"]["n_depends"], "n_pairs": 45,
            "fraction": prim["estimators"]["n_depends"] / 45.0,
            "alt_reading_defined_for_both_cells": {"n_depends": prim["estimators"]["n_depends_alt_reading"],
                                                    "fraction": prim["estimators"]["n_depends_alt_reading"] / 45.0},
            "pairs_that_depend": [r["protocol"] for r in prim["estimators"]["per_pair"] if r["depends"]],
            "min_verdicts_any_pair": min(r["n_verdicts"] for r in prim["estimators"]["per_pair"]),
        },
        "item7_fig2c_cell_to_cell_variability": {
            g: {"n_protocols": v["n_protocols"], "protocols": v["protocols"],
                "SOH_0.90": v["0.9"]["mean"], "SOH_0.875": v["0.875"]["mean"], "SOH_0.85": v["0.85"]["mean"],
                "percent_3dp": {"0.90": pct3(v["0.9"]["mean"]), "0.875": pct3(v["0.875"]["mean"]), "0.85": pct3(v["0.85"]["mean"])},
                "n_nonnan": {"0.90": v["0.9"]["n_nonnan"], "0.875": v["0.875"]["n_nonnan"], "0.85": v["0.85"]["n_nonnan"]},
                "n_one_cell_nan_counted_as_zero": {"0.90": v["0.9"]["n_protocols_with_one_cell_nan_counted_as_0"],
                                                   "0.875": v["0.875"]["n_protocols_with_one_cell_nan_counted_as_0"],
                                                   "0.85": v["0.85"]["n_protocols_with_one_cell_nan_counted_as_0"]},
                "max_over_protocols": {"0.90": v["0.9"]["max"], "0.875": v["0.875"]["max"], "0.85": v["0.85"]["max"]}}
            for g, v in fig2c.items()},
        "item7_notes": {"cells_nonmonotone_after_notebook_filter": nb_nonmono,
                        "protocols_with_one_nan_cell_on_grid": nb_one_nan},
        "item8_top5_pair_differences_at_0.90": [
            {"protocol": r[0], "cells": [r[1], r[2]], "L": [r[3], r[4]], "diff_fraction": r[5], "diff_percent_3dp": pct3(r[5])}
            for r in P[0.90]["top5"]],
        "details_per_pair": [
            {"protocol": p, "cells": [ca, cb],
             **{f"L_{s}": [prim['per_threshold'][s]['rows'][i][3], prim['per_threshold'][s]['rows'][i][4]] for s in (0.95, 0.925, 0.90, 0.875, 0.85)},
             **{f"d_{s}": prim['per_threshold'][s]['rows'][i][5] for s in (0.95, 0.925, 0.90, 0.875, 0.85)},
             "estimator_verdicts": prim["estimators"]["per_pair"][i]["verdicts"],
             "depends": prim["estimators"]["per_pair"][i]["depends"]}
            for i, (p, (ca, cb)) in enumerate([(r[0], (r[1], r[2])) for r in prim["per_threshold"][0.90]["rows"]])],
        "sensitivity_other_filters": {},
    }
    for label, res in (("runmin_from_start_keep_if_le_all_earlier", alt_b), ("notebook_adjacent_filter", alt_nb)):
        Q = res["per_threshold"]
        report["sensitivity_other_filters"][label] = {
            str(s): {"n_pairs_both": Q[s]["n_pairs_both"], "mean_pct": pct3(Q[s]["mean"]), "median_pct": pct3(Q[s]["median"]),
                     "p95_pct": pct3(Q[s]["p95"]), "max_pct": pct3(Q[s]["max"]), "max_protocol": Q[s]["max_protocol"],
                     "sigma_pct": pct3(Q[s]["sigma"]), "within_pair_share": Q[s]["within_pair_share"]}
            for s in (0.90, 0.85)}
        report["sensitivity_other_filters"][label]["estimator_dependence_fraction"] = res["estimators"]["n_depends"] / 45.0
        report["sensitivity_other_filters"][label]["top5_at_0.90"] = [(r[0], pct3(r[5])) for r in Q[0.90]["top5"]]

    def clean(o):
        if isinstance(o, dict):
            return {str(k): clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [clean(v) for v in o]
        if isinstance(o, (float, np.floating)):
            return None if math.isnan(o) else float(o)
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.bool_):
            return bool(o)
        return o

    report["definition_choices"] = [
        "Monotone filter: keep point i iff q[i] >= q[j] for all j > i (explicit definition). The paraphrase 'running minimum "
        "from the end' read literally (q[i] == min(q[i:])) would keep only points below everything after them, so it was "
        "treated as a loose paraphrase. No exact ties exist in q, so >= vs > is immaterial. 27 cells lose at least one point.",
        "Lifetime: first crossing of s*q_f[0] on the filtered curve (strictly decreasing after filtering), linear interpolation "
        "between the bracketing filtered points; identical to np.interp to 2e-13.",
        "Pairs: protocol_name groups of size 2, cells in metadata stored order (order does not affect any reported quantity).",
        "Worst-of-N: numpy.random.default_rng(PCG64), each cell draws an index uniformly with replacement from the 90-value "
        "pool; median over 10,000,000 modules per N; SE from 10 batch medians.",
        "Estimators: slopes use the first k points of the FILTERED series (every filtered series has >= 10 points, so all "
        "defined); fade at 200/400 EFC by linear interpolation of soh on the filtered series (both always within range). "
        "'At least 3 estimators defined' read as >= 3 estimators naming a worse cell; the alternative reading (>= 3 defined "
        "for both cells) gives the same count because every pair has >= 11 verdicts.",
        "Item 7: reimplemented the notebook logic and confirmed it against a literal execution of figure2.ipynb code cells "
        "1,5,6,8,9,10,11,17 (difference exactly 0).",
    ]
    with open(OUT + "recompute.json", "w") as f:
        json.dump(clean(report), f, indent=2, allow_nan=False)

    # ---------------- console summary
    print("numpy", np.__version__, "pandas", pd.__version__)
    print("pairs", prim["n_pairs"], "singles", prim["single_cell_protocols"])
    for s in (0.90, 0.85):
        q = P[s]
        print(f"s={s}: n_both={q['n_pairs_both']} mean={100*q['mean']:.6f}% median={100*q['median']:.6f}% "
              f"p95={100*q['p95']:.6f}% max={100*q['max']:.6f}% ({q['max_protocol']}) sigma={100*q['sigma']:.6f}%")
        print("   missing:", q["pairs_missing"])
    print(f"worst-of-8 median shortfall {100*w8:.6f}% (SE {100*se[8]:.5f})  worst-of-96 {100*w96:.6f}% (SE {100*se[96]:.5f})"
          f"  (seeds {SEED},{SEED+1}, {NMOD} modules)")
    print("   other seeds N8:", [f"{100*x:.4f}" for x in spread[8]], "N96:", [f"{100*x:.4f}" for x in spread[96]])
    print(f"within-pair share at 0.90: {P[0.90]['within_pair_share']:.6f}  (W={P[0.90]['W']:.3f}, T={P[0.90]['T']:.3f})")
    print(f"estimator dependence: {prim['estimators']['n_depends']}/45 = {prim['estimators']['n_depends']/45:.6f}; "
          f"alt reading {prim['estimators']['n_depends_alt_reading']}/45; min verdicts per pair "
          f"{min(r['n_verdicts'] for r in prim['estimators']['per_pair'])}")
    print("top5 at 0.90:")
    for r in P[0.90]["top5"]:
        print(f"   {r[0]:28s} {r[1]} {r[2]} L=({r[3]:.3f},{r[4]:.3f}) d={100*r[5]:.6f}%")
    print("fig2c reimplementation:")
    for g, v in fig2c.items():
        print(f"   {g}: n_prot={v['n_protocols']} 0.90={v['0.9']['mean']:.8f} 0.875={v['0.875']['mean']:.8f} "
              f"0.85={v['0.85']['mean']:.8f}  nonnan={v['0.9']['n_nonnan']},{v['0.875']['n_nonnan']},{v['0.85']['n_nonnan']}"
              f"  one-nan-as-0={v['0.9']['n_protocols_with_one_cell_nan_counted_as_0']},"
              f"{v['0.875']['n_protocols_with_one_cell_nan_counted_as_0']},{v['0.85']['n_protocols_with_one_cell_nan_counted_as_0']}")
    print("   notebook-filter non-monotone cells:", nb_nonmono)
    print("   protocols with one NaN cell on grid:", {k: (min(v), max(v), len(v)) for k, v in nb_one_nan.items()})
    for label, res in (("alt_B", alt_b), ("alt_NB", alt_nb)):
        Q = res["per_threshold"]
        print(label, {str(s): (Q[s]["n_pairs_both"], round(100*Q[s]["mean"], 3), round(100*Q[s]["median"], 3),
                              round(100*Q[s]["p95"], 3), round(100*Q[s]["max"], 3), Q[s]["max_protocol"]) for s in (0.90, 0.85)},
              "dep", res["estimators"]["n_depends"])


if __name__ == "__main__":
    main()
