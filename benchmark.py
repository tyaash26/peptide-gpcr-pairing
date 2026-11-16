#!/usr/bin/env python3
"""
Benchmark the co-evolution scores: can they rank true peptide-receptor pairs
above negatives?

This script decides nothing by eye. It reports AUROC and AUPRC with bootstrap
confidence intervals, across seeds and split schemes, plus the controls needed
to make a NULL result trustworthy -- a null is only meaningful if the measuring
instrument demonstrably works.

Controls
--------
  self-pair        every gene scored against itself. Must give ~1.0. If it
                   does not, the metric is broken and nothing else matters.
  shuffled labels  labels permuted within each fold. Must give ~0.5. If it
                   does not, the evaluation is leaking.
  presence counts  AUROC using only n_receptor or n_peptide as the score. If
                   conservation alone ranks pairs well, any apparent signal in
                   the similarity metrics may be conservation in disguise.

Usage
-----
    python benchmark.py
    python benchmark.py --metric all__bit_identity

Outputs
-------
    data/benchmark_results.csv
    data/benchmark_report.txt
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

METRICS = ["all__bit_identity", "all__pearson", "all__jaccard",
           "shared__jaccard", "shared__bit_identity", "shared__pearson"]
CONFOUNDERS = ["all__n_receptor", "all__n_peptide", "all__n_both"]


def auroc(y, s):
    """Rank-based AUROC; ties handled by average ranks. NaN scores dropped."""
    y = np.asarray(y, dtype=float)
    s = np.asarray(s, dtype=float)
    keep = ~np.isnan(s)
    y, s = y[keep], s[keep]
    n_pos, n_neg = int((y == 1).sum()), int((y == 0).sum())
    if n_pos == 0 or n_neg == 0:
        return np.nan
    order = np.argsort(s)
    ranks = np.empty(len(s), dtype=float)
    ranks[order] = np.arange(1, len(s) + 1)
    # average ranks for ties
    df = pd.DataFrame({"s": s, "r": ranks})
    ranks = df.groupby("s")["r"].transform("mean").to_numpy()
    return (ranks[y == 1].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


def auprc(y, s):
    """Average precision."""
    y = np.asarray(y, dtype=float)
    s = np.asarray(s, dtype=float)
    keep = ~np.isnan(s)
    y, s = y[keep], s[keep]
    if (y == 1).sum() == 0 or (y == 0).sum() == 0:
        return np.nan
    order = np.argsort(-s)
    y = y[order]
    tp = np.cumsum(y)
    prec = tp / np.arange(1, len(y) + 1)
    return float((prec * y).sum() / y.sum())


def boot_ci(y, s, fn, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    y = np.asarray(y)
    s = np.asarray(s, dtype=float)
    vals = []
    idx = np.arange(len(y))
    for _ in range(n):
        b = rng.choice(idx, size=len(idx), replace=True)
        v = fn(y[b], s[b])
        if not np.isnan(v):
            vals.append(v)
    if not vals:
        return np.nan, np.nan
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def self_pair_control(fp_path, say):
    """Score every gene against itself. Must be ~1.0 for identity metrics."""
    fp = pd.read_csv(fp_path, index_col=0)
    mat = fp.drop(columns=[c for c in ("role",) if c in fp.columns])
    mat = mat.astype(float)
    nonzero = mat[mat.sum(axis=1) > 0]
    bits, jac = [], []
    for i in nonzero.index:
        v = nonzero.loc[i].to_numpy()
        bits.append(float((v == v).mean()))
        both = float(((v == 1) & (v == 1)).sum())
        either = float(((v == 1) | (v == 1)).sum())
        jac.append(both / either if either else np.nan)
    say(f"  self-pair bit_identity: {np.mean(bits):.4f} "
        f"(expected 1.0000)")
    say(f"  self-pair jaccard:      {np.nanmean(jac):.4f} "
        f"(expected 1.0000)")
    ok = abs(np.mean(bits) - 1.0) < 1e-9 and abs(np.nanmean(jac) - 1.0) < 1e-9
    say(f"  instrument check: {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", default="data")
    ap.add_argument("--fingerprints", default="data/fingerprints_all.csv")
    ap.add_argument("--out", default="data")
    ap.add_argument("--boot", type=int, default=2000)
    a = ap.parse_args()

    log = []

    def say(m):
        print(m)
        log.append(m)

    say("CONTROL 1: self-pairs (does the metric detect identity at all?)")
    self_pair_control(a.fingerprints, say)

    files = sorted(glob.glob(os.path.join(a.scores, "scores_*_seed*.csv")))
    if not files:
        sys.exit(f"no score files under {a.scores}")
    say(f"\nscore files: {len(files)}")

    rows = []
    for path in files:
        m = re.match(r"scores_(.+)_(hard|random|matched)_seed(\d+)\.csv",
                     os.path.basename(path))
        if not m:
            continue
        scheme, negs, seed = m.groups()
        d = pd.read_csv(path)
        if d.empty or d["label"].nunique() < 2:
            continue
        y = d["label"].to_numpy()

        present = [c for c in METRICS + CONFOUNDERS if c in d.columns]
        for metric in present:
            s = d[metric].to_numpy(dtype=float)
            ro, pr = auroc(y, s), auprc(y, s)
            lo, hi = boot_ci(y, s, auroc, n=a.boot, seed=int(seed))
            rows.append({"scheme": scheme, "negatives": negs,
                         "seed": int(seed), "metric": metric,
                         "n": len(d), "n_pos": int((y == 1).sum()),
                         "auroc": ro, "auroc_lo": lo, "auroc_hi": hi,
                         "auprc": pr, "base_rate": float((y == 1).mean()),
                         "kind": "confounder" if metric in CONFOUNDERS
                                 else "similarity"})

        # CONTROL 2: shuffled labels on the primary metric
        rng = np.random.default_rng(int(seed) + 99)
        ysh = rng.permutation(y)
        if "all__bit_identity" in d.columns:
            s = d["all__bit_identity"].to_numpy(dtype=float)
            rows.append({"scheme": scheme, "negatives": negs,
                         "seed": int(seed), "metric": "SHUFFLED_LABELS",
                         "n": len(d), "n_pos": int((ysh == 1).sum()),
                         "auroc": auroc(ysh, s), "auroc_lo": np.nan,
                         "auroc_hi": np.nan, "auprc": auprc(ysh, s),
                         "base_rate": float((ysh == 1).mean()),
                         "kind": "control"})

    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(a.out, "benchmark_results.csv"), index=False)

    # --- summary across seeds -------------------------------------------
    say("\nAUROC by scheme and metric (mean across seeds, with per-seed range)")
    for (scheme, negs), g in res.groupby(["scheme", "negatives"]):
        say(f"\n{scheme} / {negs} negatives  "
            f"(n={int(g['n'].mean())}, base rate {g['base_rate'].mean():.2f})")
        for kind in ("similarity", "confounder", "control"):
            sub = g[g["kind"] == kind]
            if sub.empty:
                continue
            for metric, gm in sub.groupby("metric"):
                mu = gm["auroc"].mean()
                lo, hi = gm["auroc"].min(), gm["auroc"].max()
                ap_ = gm["auprc"].mean()
                tag = {"confounder": "  [confounder]",
                       "control": "  [control]"}.get(kind, "")
                say(f"   {metric:<24} AUROC {mu:.3f} "
                    f"(seeds {lo:.3f}-{hi:.3f})  AUPRC {ap_:.3f}{tag}")

    say("")
    say("How to read this. A similarity metric shows signal only if its AUROC "
        "is clearly above 0.5 with a bootstrap interval excluding 0.5, while "
        "SHUFFLED_LABELS sits at 0.5. If a confounder (n_receptor, n_peptide) "
        "scores as well as the similarity metrics, the metrics are tracking "
        "how conserved the genes are, not whether they are partners. The "
        "hard-negative columns are the primary result; random negatives are "
        "shown only for contrast with published comparisons.")

    print(f"\nwrote {a.out}/benchmark_results.csv")
    with open(os.path.join(a.out, "benchmark_report.txt"), "w") as f:
        f.write("\n".join(log) + "\n")
    print("wrote", os.path.join(a.out, "benchmark_report.txt"))


if __name__ == "__main__":
    main()
