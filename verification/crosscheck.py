"""Run:  python3 -I verification/crosscheck.py <path to dynamic_cycling_Nature_Energy_2024/data>
It executes the computational code cells of the authors' figure2.ipynb.

Cross-checks: (a) literal execution of figure2.ipynb computational code cells vs my reimplementation,
(b) lifetime() vs np.interp on primary-filtered series, (c) high-precision worst-of-N medians."""
import json
import sys

import numpy as np
import pandas as pd

import os  # noqa: E402
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import recompute as R  # noqa: E402

NB = os.path.join(os.path.dirname(R.DATA.rstrip("/")), "figure2.ipynb")

meta, diag = R.verify_and_load()
series = R.cell_series(diag)

# ---------------------------------------------------------------- (a) literal notebook execution
with open(NB, "r", encoding="utf-8") as f:
    nb = json.load(f)
wanted = [1, 5, 6, 8, 9, 10, 11, 17]
ns = {"np": np, "pd": pd}
for i in wanted:
    c = nb["cells"][i]
    assert c["cell_type"] == "code"
    src = "".join(c["source"]) if isinstance(c["source"], list) else c["source"]
    src = src.replace("'./data/", "'" + R.DATA)
    exec(compile(src, f"figure2_cell{i}", "exec"), ns)

ivm = ns["intr_var_mean"]
mine, nonmono, _ = R.fig2c_reimplementation(meta, series)
print("literal notebook exec, intr_var_mean (C/10, C/5, C/2):")
maxabs = 0.0
for soh in (0.90, 0.875, 0.85):
    lit = ivm[soh]
    me = [mine[g][f"{soh}"]["mean"] for g in ("C/10", "C/5", "C/2")]
    maxabs = max(maxabs, max(abs(a - b) for a, b in zip(lit, me)))
    print(f"  SOH {soh}: notebook {[f'{x:.10f}' for x in lit]}  mine {[f'{x:.10f}' for x in me]}")
print("  max |notebook - mine| =", maxabs)
for cr in ns["crates"]:
    pr = ns["data"][ns["data"]["cell_name"].isin(ns["buckets3"][cr])]["protocol_name"].unique()
    pr = pr[(pr != "CC_(no_storage)_Co2") & (pr != "Synthetic_2c_Co5")]
    print(f"  {cr}: {len(pr)} protocols; same set as mine: {sorted(pr) == sorted(mine['C/' + cr[2:]]['protocols'])}")
# full-grid equality of the per-cell EFC table
tab = ns["allSOHs_pd"]
SOHs = np.linspace(1, 0.85, 31)
worst = 0.0
for cell in tab.index:
    q, e = series[cell]
    qf, ef = R.filter_notebook(q, e)
    v = np.round(np.interp(SOHs * qf[0], np.flip(qf), np.flip(ef), left=np.nan, right=np.nan), 1)
    v[0] = 0
    a = tab.loc[cell].to_numpy()
    both = ~(np.isnan(a) & np.isnan(v))
    assert np.array_equal(np.isnan(a), np.isnan(v)), cell
    worst = max(worst, float(np.nanmax(np.abs(a[both] - v[both])) if both.any() else 0.0))
print("  allSOHs table identical to mine (max abs diff):", worst)

# what np.interp returns on the two non-monotone notebook series vs a first-crossing reading
for cell in nonmono:
    q, e = series[cell]
    qf, ef = R.filter_notebook(q, e)
    for soh in (0.90, 0.875, 0.85):
        j = int(np.where(np.round(SOHs, 3) == soh)[0][0])
        print(f"  {cell} SOH {soh}: notebook grid value {tab.loc[cell].iloc[j]}  first-crossing on same series "
              f"{R.lifetime(qf, ef, soh):.3f}")

# sensitivity: Fig 2c with the primary (strict) filter in place of the notebook filter
orig = R.filter_notebook
R.filter_notebook = R.filter_primary
alt, _, _ = R.fig2c_reimplementation(meta, series)
R.filter_notebook = orig
print("  Fig2c if the primary filter replaced the notebook filter:",
      {g: [round(alt[g][f'{s}']['mean'], 6) for s in (0.90, 0.875, 0.85)] for g in alt})

# ---------------------------------------------------------------- (b) lifetime vs np.interp
dev = 0.0
for cell, (q, e) in series.items():
    qf, ef = R.filter_primary(q, e)
    assert np.all(np.diff(qf) < 0) and np.all(np.diff(ef) > 0)
    for s in R.LIFE_S:
        a = R.lifetime(qf, ef, s)
        b = float(np.interp(s * qf[0], qf[::-1], ef[::-1], left=np.nan, right=np.nan))
        if np.isnan(a) or np.isnan(b):
            assert np.isnan(a) and np.isnan(b), (cell, s, a, b)
        else:
            dev = max(dev, abs(a - b))
print("(b) max |lifetime - np.interp| over all cells and thresholds:", dev)

# ---------------------------------------------------------------- (c) high-precision worst-of-N
prim, _ = R.analyse(meta, series, R.filter_primary, with_estimators=False)
d90 = prim["per_threshold"][0.90]["d"]
for N, seed in ((8, 777), (96, 778)):
    med, _ = R.worst_of_n_median(d90, N, 10_000_000, seed, chunk=100_000)
    print(f"(c) N={N}: median shortfall over 10,000,000 modules (seed {seed}) = {100*med:.5f}%")
