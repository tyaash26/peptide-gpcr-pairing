#!/usr/bin/env python3
"""
Permutation null for the benchmark, replacing the single shuffled-label check.

Why
---
benchmark.py shuffles labels once per seed as a control, and on the
class-matched splits that single shuffle returned AUROC 0.458 (range
0.451-0.473) rather than ~0.5. At n = 571 a single shuffle is a noisy estimate
of the null, and it is off by roughly the same magnitude as the apparent
signal (0.529-0.559). Comparing an observed AUROC to a theoretical 0.5 is
therefore not safe here.

This builds the null empirically: permute the labels many times, record the
AUROC each time, and report where the observed value falls in that
distribution. The p-value is the fraction of permutations at least as extreme.

Permutation respects the fold structure: labels are shuffled WITHIN each fold,
so any fold-level imbalance in base rate is preserved in the null rather than
contributing to it.

Also reports per-family AUROC, because the recurring confound in this project
is a metric separating receptor families rather than detecting pairs.

Usage
-----
    python permutation_test.py --scores data/class_scored
    python permutation_test.py --scores data/class_scored --n-perm 20000

Outputs
-------
    <scores>/permutation_results.csv
    <scores>/permutation_report.txt
"""

import argparse
import glob
import importlib.util
import os
import re
import sys

import numpy as np
import pandas as pd

METRICS = ["all__bit_identity", "all__jaccard", "all__pearson",
           "shared__jaccard", "shared__pearson",
           "all__n_receptor", "all__n_peptide"]


def load_bench(path="benchmark.py"):
    if not os.path.exists(path):
        sys.exit(f"{path} not found; run from the project root")
    spec = importlib.util.spec_from_file_location("bench", path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def permute_within_folds(y, folds, rng):
    out = y.copy()
    for f in np.unique(folds):
        idx = np.where(folds == f)[0]
        out[idx] = rng.permutation(y[idx])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", default="data/class_scored")
    ap.add_argument("--pairs", default="data/pairs_resolved.csv")
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    b = load_bench()
    log = []

    def say(m=""):
        print(m)
        log.append(m)

    files = sorted(glob.glob(os.path.join(a.scores, "scores_*_seed*.csv")))
    if not files:
        sys.exit(f"no score files under {a.scores}")

    pairs = pd.read_csv(a.pairs, dtype=str)
    fam_of = dict(zip(pairs["receptor_uniprot"], pairs["receptor_family"]))

    rows = []
    for path in files:
        m = re.match(r"scores_(.+)_(\w+)_seed(\d+)\.csv", os.path.basename(path))
        if not m:
            continue
        scheme, negs, seed = m.groups()
        d = pd.read_csv(path)
        if d.empty or d["label"].nunique() < 2:
            continue
        y = d["label"].to_numpy()
        folds = d["fold"].to_numpy() if "fold" in d.columns else np.zeros(len(d))

        say(f"\n{os.path.basename(path)}  n={len(d)}, "
            f"base rate {y.mean():.3f}, {a.n_perm:,} permutations")

        for metric in METRICS:
            if metric not in d.columns:
                continue
            s = d[metric].to_numpy(dtype=float)
            obs = b.auroc(y, s)

            rng = np.random.default_rng(a.seed)
            null = np.empty(a.n_perm)
            for i in range(a.n_perm):
                null[i] = b.auroc(permute_within_folds(y, folds, rng), s)
            null = null[~np.isnan(null)]

            centre = null.mean()
            lo, hi = np.percentile(null, [2.5, 97.5])
            # two-sided p: how often the null is at least as far from its own
            # centre as the observed value is
            p = (np.abs(null - centre) >= abs(obs - centre)).mean()
            p = max(p, 1.0 / len(null))

            tag = "  [confounder]" if metric.endswith(
                ("n_receptor", "n_peptide")) else ""
            flag = "  <-- exceeds null" if (obs > hi or obs < lo) else ""
            say(f"  {metric:<22} obs {obs:.3f}  null centre {centre:.3f} "
                f"[{lo:.3f}–{hi:.3f}]  p={p:.4f}{tag}{flag}")
            rows.append({"file": os.path.basename(path), "seed": int(seed),
                         "metric": metric, "observed": obs,
                         "null_centre": centre, "null_lo": lo, "null_hi": hi,
                         "p_value": p, "n": len(d)})

        # per-family check for the strongest similarity metric
        if ("all__bit_identity" in d.columns and "receptor_uniprot" in d.columns
                and seed == "0"):
            d = d.copy()
            d["fam"] = d["receptor_uniprot"].map(fam_of)
            say("  per-family AUROC (all__bit_identity, families with >= 40 rows)")
            any_fam = False
            for fam, g in d.groupby("fam"):
                if len(g) < 40 or g["label"].nunique() < 2:
                    continue
                any_fam = True
                say(f"    {str(fam)[:36]:<36} n={len(g):>4} "
                    f"pos={int(g.label.sum()):>3} "
                    f"AUROC={b.auroc(g.label.to_numpy(), g['all__bit_identity'].to_numpy()):.3f}")
            if not any_fam:
                say("    no family has enough rows to evaluate")

    if not rows:
        sys.exit("nothing evaluated")

    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(a.scores, "permutation_results.csv"), index=False)

    say()
    say("Reading this. The null centre is where a no-signal result actually "
        "sits for this data, which need not be 0.500. An observed AUROC only "
        "counts as signal if it falls outside the null interval and the "
        "confounders do not. If the per-family AUROCs sit at chance while the "
        "pooled value does not, the metric is separating receptor families "
        "rather than identifying pairs.")

    with open(os.path.join(a.scores, "permutation_report.txt"), "w") as f:
        f.write("\n".join(log) + "\n")
    print(f"\nwrote {a.scores}/permutation_results.csv and "
          f"{a.scores}/permutation_report.txt")


if __name__ == "__main__":
    main()
