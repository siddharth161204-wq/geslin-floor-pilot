#!/usr/bin/env python3
"""
Synthetic check of lifetime_floor.py, run before the real analysis.

It keeps the real experimental design (the 92 cells and 47 protocols of data/metadata.pkl) and
replaces every capacity trajectory with a synthetic one whose per-cell lifetime noise is known,
then runs lifetime_floor.py on that folder and checks that the script recovers what was put in:

  1. the per-cell noise sigma (true value 3 %, plus a small known contribution from capacity
     measurement noise), within the sampling band that 45 pairs allow;
  2. the normal-model worst-of-N medians equal the normal-theory multiples of the fitted sigma
     (1.365 sigma for N = 8, 2.445 sigma for N = 96) to Monte Carlo accuracy;
  3. the nonparametric pool has the variance of a single cell (mean of delta squared over 2);
  4. the within-pair variance share is close to its known value;
  5. the oriented pair ratio crosses 1 for exactly the pairs whose worse cell depends on the
     estimator;
  6. the notebook's single-pass filter and the primary filter disagree only where the notebook's
     output is non-monotone, and the reconciliation block runs.

No real capacity value is used. Usage:
  python3 tests/synthetic_check.py --data /path/to/dynamic_cycling_Nature_Energy_2024/data
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, os.pardir, "lifetime_floor.py")
sys.path.insert(0, os.path.join(HERE, os.pardir))
import lifetime_floor as lf  # noqa: E402

SIGMA_TRUE = 0.03       # per-cell relative lifetime noise put into the synthetic data
NOISE_Q = 0.0008        # capacity measurement noise (fraction of initial capacity)
BUMP_P, BUMP = 0.06, 0.004


def g(u):
    """Normalised fade shape: SOH = 1 - 0.10 g(EFC / L90), g(1) = 1."""
    return 0.35 * np.sqrt(u) + 0.65 * u


def synthesise(meta, rng):
    rows, eol = [], {}
    protocols = meta["protocol_name"].unique()
    l90_p = dict(zip(protocols, rng.uniform(350.0, 1300.0, len(protocols))))
    short = {"CC_(5%_storage)_Co10", "CC_(50%_storage)_Co10"}     # stopped before 85 % SOH
    for cell, prot in zip(meta["cell_name"], meta["protocol_name"]):
        l90 = l90_p[prot] * (1.0 + rng.normal(0.0, SIGMA_TRUE))
        cycles = [0, 25, 50, 75, 100]
        # 85 % SOH is reached at about 1.62 L90; the short protocols stop between 90 and 85 %.
        efc_stop = (1.35 * l90_p[prot] + 100.0) if prot in short else 2.4 * l90_p[prot]
        while True:
            k = cycles[-1] + 100
            if 0.6 + 0.95 * k + 1.3 * len(cycles) > efc_stop:
                break
            cycles.append(k)
        efc = np.array([0.6 + 0.95 * k + 1.3 * i for i, k in enumerate(cycles)])
        soh = 1.0 - 0.10 * g(efc / l90)
        q = soh + rng.normal(0.0, NOISE_Q, len(efc))
        q[1:] += (rng.random(len(efc) - 1) < BUMP_P) * BUMP
        if cell == "cell_077":
            # At the 90 % SOH crossing: the capacity dips below the threshold for two diagnostics and
            # recovers above it. The notebook's single pass keeps the first dipped point, so its curve
            # crosses the threshold twice; the monotone filter drops both dipped points.
            t = 0.9 * q[0]
            j = min(max(int(np.argmax(q < t)), 1), len(q) - 4)
            q[j - 1] = max(q[j - 1], t + 0.003)
            q[j], q[j + 1], q[j + 2] = t - 0.002, t - 0.004, t + 0.001
            q[j + 3] = min(q[j + 3], t - 0.006)
        if cell == "cell_085":
            # mid-record, away from any threshold: non-monotone under the notebook, same lifetimes
            j = len(q) // 2
            q[j] = q[j - 1] - 0.006
            q[j + 1] = q[j - 1] + 0.001
        for k, e, v in zip(cycles, efc, q):
            rows.append((cell, float(k), float(v), float(e)))
        eol[cell] = l90 + 0.6      # what a perfect paper-side lifetime at 90 % would be, roughly
    dfa = pd.DataFrame(rows, columns=["Cell", "Cycle", lf.CAPA, lf.EFC]).set_index(["Cell", "Cycle"])
    em = pd.DataFrame({"EFCs (with Diagnostic)": pd.Series(eol)})
    return dfa, em


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="the real data folder (only metadata.pkl is read)")
    ap.add_argument("--keep", default=None, help="keep the synthetic folder and results here")
    args = ap.parse_args()
    rng = np.random.default_rng(12345)
    meta = lf.read_pickle(os.path.join(args.data, "metadata.pkl"))
    work = args.keep or tempfile.mkdtemp(prefix="lifetime_floor_synthetic_")
    syn = os.path.join(work, "data")
    out = os.path.join(work, "results")
    os.makedirs(syn, exist_ok=True)
    shutil.copy(os.path.join(args.data, "metadata.pkl"), os.path.join(syn, "metadata.pkl"))
    dfa, em = synthesise(meta, rng)
    dfa.to_pickle(os.path.join(syn, "diagnostic_features_all.pkl"))
    em.to_pickle(os.path.join(syn, "eol_metrics_soh0.9.pkl"))
    subprocess.run([sys.executable, SCRIPT, "--data", syn, "--out", out, "--allow-unverified-inputs"],
                   check=True, stdout=subprocess.DEVNULL)
    r = json.load(open(os.path.join(out, "results.json")))

    checks = []

    def check(name, ok, detail):
        checks.append((name, bool(ok), detail))

    # The loader must refuse a pickle that would run code, whatever its module prefix.
    class _Evil:
        def __init__(self, fn, arg):
            self.fn, self.arg = fn, arg

        def __reduce__(self):
            return (self.fn, (self.arg,))

    import pickle
    import numpy.testing._private.utils as npu
    refused = []
    for fn, arg in ((os.system, "echo SHOULD-NOT-RUN"), (npu.runstring, "raise SystemExit(99)")):
        path = os.path.join(work, "evil.pkl")
        with open(path, "wb") as f:
            f.write(pickle.dumps(_Evil(fn, arg)))
        try:
            lf.read_pickle(path)
            refused.append(False)
        except pickle.UnpicklingError:
            refused.append(True)
    check("loader refuses code-running globals (os.system, numpy runstring)", all(refused), str(refused))

    a2 = r["A2_cell_noise_sigma"]["soh0.900"]
    sig = a2["sigma_from_mean"]
    # Expected: true lifetime noise plus capacity noise converted to lifetime through the slope at
    # 90 % SOH (0.10 g'(1) per unit relative life), shared between interpolation neighbours.
    meas = NOISE_Q / (0.10 * (0.35 / 2 + 0.65))
    check("design: 45 pairs at 90 %, singletons cell_017 and cell_084",
          r["sanity"]["n_pairs_at_soh90"] == 45 and sorted(s["cell"] for s in r["singletons"]) == ["cell_017", "cell_084"],
          f"{r['sanity']['n_pairs_at_soh90']} pairs")
    check("sigma recovered within the band 45 pairs allow", 0.022 <= sig <= 0.042,
          f"sigma_from_mean {sig:.4f} against true {SIGMA_TRUE} (plus up to about {meas:.4f} from capacity noise)")
    w = r["A3_selection_only_worst_of_N"]["soh0.900"]["worst_of_N_shortfall"]
    for n, mult in (("8", 1.3647), ("96", 2.4452)):
        got = w[n]["parametric_median"] / sig
        check(f"normal-model worst-of-{n} median equals {mult} sigma", abs(got / mult - 1) < 0.02, f"{got:.4f} sigma")
    pt = pd.read_csv(os.path.join(out, "pair_differences.csv"))
    d = pt.loc[np.isclose(pt["threshold"], 0.90), "delta_rel"].to_numpy()
    pool = lf.deviation_pool(d)
    check("pool variance equals mean(delta^2)/2", abs(np.mean(pool ** 2) / (np.mean(d ** 2) / 2) - 1) < 1e-12,
          f"{np.mean(pool ** 2):.3e} vs {np.mean(d ** 2) / 2:.3e}")
    share = r["A5_variance_split"]["soh0.900"]["within_pair_variance_share"]
    check("within-pair share small and positive", 0.0 < share < 0.05, f"{share:.4f}")
    rr = pd.read_csv(os.path.join(out, "pair_ratio_ranges.csv"))
    ok = rr[rr["n_estimators"] >= 3]
    check("oriented ratio crosses 1 exactly where the worse cell flips",
          bool(((ok["oriented_ratio_min"] < 1.0) == ok["worse_cell_depends_on_estimator"]).all()),
          f"{int((ok['oriented_ratio_min'] < 1).sum())} crossing, {int(ok['worse_cell_depends_on_estimator'].sum())} flipping")
    rec = r["R_reconciliation"]["paper_fig2c"]
    check("reconciliation groups are 13, 9 and 6 protocols",
          [rec["groups"][k]["n_protocols"] for k in ("C/10", "C/5", "C/2")] == [13, 9, 6],
          str([rec["groups"][k]["n_protocols"] for k in ("C/10", "C/5", "C/2")]))
    nonmono = []
    for cell in meta["cell_name"]:
        sub = dfa.loc[cell]
        e, q = lf.notebook_filter(np.asarray(sub[lf.EFC]), np.asarray(sub[lf.CAPA]))
        if np.any(np.diff(q) > 0):
            nonmono.append(cell)
    differing = set()
    for blk in rec["primary_versus_notebook_pair_differences"]:
        differing |= set(blk["pairs_differing_by_more_than_0.001"])
    unexplained = []
    for prot in sorted(differing):
        row = pt[pt["protocol"] == prot].iloc[0]
        if row["cell_a"] not in nonmono and row["cell_b"] not in nonmono:
            unexplained.append(prot)
    prot77 = meta.loc[meta["cell_name"] == "cell_077", "protocol_name"].iloc[0]
    sub = dfa.loc["cell_077"]
    e77, q77 = np.asarray(sub[lf.EFC], float), np.asarray(sub[lf.CAPA], float)
    # The notebook evaluates the whole SOH grid in one vectorised np.interp call; on a non-monotone
    # curve the interval numpy's search lands in can depend on the previous grid point, so the
    # comparison is made the notebook's way.
    grid = np.linspace(1, 0.85, 31)
    pick = [20, 25, 30]                                    # 0.90, 0.875, 0.85
    em_, qm_ = lf.monotone_filter(e77, q77)
    en_, qn_ = lf.notebook_filter(e77, q77)
    l_mono = np.interp(grid * qm_[0], np.flip(qm_), np.flip(em_), left=np.nan, right=np.nan)[pick]
    l_nb = np.interp(grid * qn_[0], np.flip(qn_), np.flip(en_), left=np.nan, right=np.nan)[pick]
    gap = np.nanmax(np.abs(l_mono - l_nb) / l_mono)
    check("the injected dip gives a different lifetime under the two filters (test is not vacuous)", gap > 0.002,
          f"cell_077 at 90/87.5/85 %: monotone {np.round(l_mono, 1).tolist()}, notebook {np.round(l_nb, 1).tolist()} EFC")
    check("filters disagree only in pairs with a non-monotone notebook curve, and the injected pair is flagged",
          not unexplained and prot77 in differing and sorted(nonmono) == ["cell_077", "cell_085"],
          f"non-monotone: {sorted(nonmono)}; pairs differing: {sorted(differing)}; unexplained: {unexplained}")
    check("figure written", os.path.exists(os.path.join(out, "lifetime_floor.png")), "lifetime_floor.png")

    width = max(len(c[0]) for c in checks)
    for name, ok_, detail in checks:
        print(f"{'PASS' if ok_ else 'FAIL'}  {name.ljust(width)}  {detail}")
    if not args.keep:
        shutil.rmtree(work)
    sys.exit(0 if all(c[1] for c in checks) else 1)


if __name__ == "__main__":
    main()
