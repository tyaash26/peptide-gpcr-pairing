#!/usr/bin/env python3
"""
Degree-matched hard negatives, to remove the leak in family-matched sampling.

The problem
-----------
make_splits.py draws hard negatives as same-family receptors a peptide does not
bind. A receptor with many known partners is therefore more likely to be the
positive for any peptide in its family, so receptor degree predicts the label
by construction: AUROC 0.690, above every similarity metric and above receptor
conservation.

The fix
-------
For each positive (peptide, receptor), draw its negative from same-family
non-partner receptors whose DEGREE is as close as possible to the positive
receptor's degree. If degree is matched between the classes, degree cannot
discriminate, and whatever AUROC the similarity metrics then achieve is not
inherited from the sampling.

Cost: the hard-negative pool is small (1,041 candidates), so matching cannot be
exact everywhere. The report states the achieved degree balance and how many
positives got no usable match, rather than hiding either.

Usage
-----
    python make_splits_balanced.py
    python make_splits_balanced.py --tolerance 2 --folds 5

Outputs (in --out, default data/splits_balanced/):
    splits_cold_receptor_matched_seed<N>.csv
    splits_balanced_report.txt
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

RECEPTOR = "receptor_gtopdb_id"
PEPTIDE = "peptide_gtopdb_id"
FAMILY = "receptor_family"


def auroc(y, s):
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
    ranks = pd.DataFrame({"s": s, "r": ranks}).groupby("s")["r"] \
             .transform("mean").to_numpy()
    return (ranks[y == 1].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default="data/pairs_resolved.csv")
    ap.add_argument("--out", default="data/splits_balanced")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--tolerance", type=int, default=1,
                    help="max |degree difference| accepted as a match")
    ap.add_argument("--per-positive", type=int, default=1,
                    help="how many degree-matched negatives to draw per "
                         "positive; more rows means narrower permutation nulls")
    ap.add_argument("--scope", choices=["family", "class"], default="family",
                    help="negative pool: same receptor family (strict control, "
                         "but unobtainable for ~52%% of pairs) or any peptide "
                         "receptor in the set (larger pool, weaker control)")
    a = ap.parse_args()

    log = []

    def say(m=""):
        print(m)
        log.append(m)

    pairs = pd.read_csv(a.pairs, dtype=str).dropna(
        subset=[RECEPTOR, PEPTIDE]).drop_duplicates(subset=[PEPTIDE, RECEPTOR])
    say(f"positives: {len(pairs):,} pairs, "
        f"{pairs[RECEPTOR].nunique()} receptors, "
        f"{pairs[PEPTIDE].nunique()} peptides")

    degree = pairs.groupby(RECEPTOR)[PEPTIDE].nunique().to_dict()
    fam = dict(zip(pairs[RECEPTOR], pairs[FAMILY]))
    by_family = {}
    for r, f in fam.items():
        by_family.setdefault(f, []).append(r)
    partners = pairs.groupby(PEPTIDE)[RECEPTOR].apply(set).to_dict()

    all_receptors = sorted(degree)
    say(f"receptor degree: median {int(np.median(list(degree.values())))}, "
        f"max {max(degree.values())}")
    say(f"negative pool scope: {a.scope}"
        + (" (same receptor family only)" if a.scope == "family"
           else f" (any of {len(all_receptors)} peptide receptors)"))
    if a.scope == "class":
        say("  NOTE: a class-wide pool is a weaker control than same-family. "
            "Negatives can differ from the positive in receptor family, so a "
            "metric that separates families will look like signal. Compare "
            "against the family-scoped run rather than reading this alone.")

    os.makedirs(a.out, exist_ok=True)

    for seed in a.seeds:
        rng = np.random.default_rng(seed)
        rows, unmatched, short, used = [], 0, 0, {}
        # one negative per positive, degree-matched within family
        for _, pos in pairs.iterrows():
            p, r = pos[PEPTIDE], pos[RECEPTOR]
            d_target = degree[r]
            if a.scope == "family":
                candidates = by_family.get(fam.get(r), [])
            else:
                candidates = all_receptors
            pool = [x for x in candidates
                    if x not in partners.get(p, set()) and x != r]
            if not pool:
                unmatched += 1
                continue
            # prefer unused candidates so negatives are not dominated by a
            # handful of low-degree receptors
            # candidates within tolerance, nearest degree first
            diffs = sorted((abs(degree[x] - d_target), x) for x in pool)
            ok = [x for d, x in diffs if d <= a.tolerance]
            if not ok:
                unmatched += 1
                continue
            take = min(a.per_positive, len(ok))
            if take < len(ok):
                # among the acceptable ones, prefer the closest degrees but
                # break ties at random so negatives are not always the same
                cutoff = sorted(abs(degree[x] - d_target) for x in ok)[take - 1]
                tier = [x for x in ok if abs(degree[x] - d_target) <= cutoff]
                picks = list(rng.choice(tier, size=take, replace=False))
            else:
                picks = ok
            if take < a.per_positive:
                short += 1
            for pick in picks:
                pick = str(pick)
                used[(p, pick)] = used.get((p, pick), 0) + 1
                rows.append({PEPTIDE: p, RECEPTOR: pick,
                             FAMILY: fam.get(pick), "label": 0})

        neg = pd.DataFrame(rows)
        pos_df = pairs[[PEPTIDE, RECEPTOR, FAMILY]].copy()
        pos_df["label"] = 1
        allp = pd.concat([pos_df, neg], ignore_index=True)
        allp["degree"] = allp[RECEPTOR].map(degree).fillna(0)

        # cold-receptor folds
        uniq = allp[RECEPTOR].unique()
        order = rng.permutation(len(uniq))
        fold_of = {g: int(order[i] % a.folds) for i, g in enumerate(uniq)}
        allp["fold"] = allp[RECEPTOR].map(fold_of)

        # leakage assertion
        for k in range(a.folds):
            te = set(allp.loc[allp.fold == k, RECEPTOR])
            tr = set(allp.loc[allp.fold != k, RECEPTOR])
            if te & tr:
                sys.exit(f"LEAKAGE seed {seed} fold {k}")

        if seed == a.seeds[0]:
            say()
            say(f"seed {seed}: {len(pos_df):,} positives, {len(neg):,} "
                f"negatives, {unmatched} positives unmatched "
                f"(tolerance {a.tolerance}); {short} positives got fewer "
                f"than {a.per_positive} negatives")
            dp = allp.loc[allp.label == 1, "degree"]
            dn = allp.loc[allp.label == 0, "degree"]
            say(f"  mean receptor degree: positives {dp.mean():.2f}, "
                f"negatives {dn.mean():.2f}")
            say(f"  median:               positives {dp.median():.0f}, "
                f"negatives {dn.median():.0f}")
            ro = auroc(allp.label.to_numpy(), allp.degree.to_numpy())
            say(f"  AUROC of receptor degree: {ro:.3f}  "
                f"(was 0.690 unmatched; 0.5 means the leak is closed)")
            say(f"  distinct receptors used as negatives: "
                f"{neg[RECEPTOR].nunique()}")
            say("  leakage check passed (cold_receptor)")

        allp.to_csv(os.path.join(
            a.out, f"splits_cold_receptor_matched_seed{seed}.csv"), index=False)

    say()
    say("Next: score these with score_coevolution.py --splits "
        f"{a.out}, then benchmark.py. If the similarity metrics still sit at "
        "chance with degree neutralised, the null result holds under a "
        "corrected design.")

    with open(os.path.join(a.out, "splits_balanced_report.txt"), "w") as f:
        f.write("\n".join(log) + "\n")
    print(f"\nwrote splits and report to {a.out}/")


if __name__ == "__main__":
    main()
