#!/usr/bin/env python3
"""
Build evaluation splits for the peptide-GPCR pairing benchmark.

No features, no model. This script only decides what is a positive, what is a
negative, and which rows may be seen during training. Every later number
depends on these choices, so they are explicit and logged.

Negatives
---------
A non-pair is NOT a verified negative: it is an untested pair. Some are real
pairings nobody has assayed. This caps achievable AUROC and is a stated
limitation, not a bug.

  hard    (peptide, receptor) where the receptor shares a receptor_family with
          a known partner of that peptide, but is not itself a known partner.
          This is the benchmark that matters.
  random  (peptide, receptor) drawn from all non-pairs. Weak baseline, kept
          only to quantify how much easier it is.

Splits
------
  cold_receptor  whole receptors held out with all their edges. This is the
                 orphan condition: the model has never seen the receptor.
                 Peptides may appear in train and test, which is correct --
                 when scoring an orphan you do know the peptides.
  cold_both      additionally holds out the peptides, a stricter condition.
  random_edge    random split over pairs. Inflated; reported for contrast.

Usage
-----
    python make_splits.py                       # 5 folds, 3 seeds, 1:1
    python make_splits.py --neg-ratio 5 --folds 5

Outputs to --out (default data/splits/):
    splits_<scheme>_<negatives>_seed<N>.csv   pair, label, fold
    splits_report.txt
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

RECEPTOR = "receptor_gtopdb_id"
PEPTIDE = "peptide_gtopdb_id"
FAMILY = "receptor_family"


def load_pairs(path):
    df = pd.read_csv(path, dtype=str)
    missing = [c for c in (RECEPTOR, PEPTIDE, FAMILY) if c not in df.columns]
    if missing:
        sys.exit(f"{path}: missing columns {missing}")
    df = df.dropna(subset=[RECEPTOR, PEPTIDE]).drop_duplicates(
        subset=[PEPTIDE, RECEPTOR])
    return df


def build_negatives(pairs, kind, ratio, rng):
    """Return a DataFrame of negative (peptide, receptor) pairs."""
    receptors = pairs[RECEPTOR].unique()
    peptides = pairs[PEPTIDE].unique()
    known = set(zip(pairs[PEPTIDE], pairs[RECEPTOR]))

    fam = dict(zip(pairs[RECEPTOR], pairs[FAMILY]))
    by_family = {}
    for r in receptors:
        by_family.setdefault(fam.get(r), []).append(r)

    # receptors each peptide is known to bind
    partners = pairs.groupby(PEPTIDE)[RECEPTOR].apply(list).to_dict()

    candidates = []
    if kind == "hard":
        for p, rs in partners.items():
            fams = {fam.get(r) for r in rs}
            pool = {r for f in fams for r in by_family.get(f, [])}
            pool -= set(rs)
            for r in pool:
                if (p, r) not in known:
                    candidates.append((p, r))
    elif kind == "random":
        for p in peptides:
            for r in receptors:
                if (p, r) not in known:
                    candidates.append((p, r))
    else:
        sys.exit(f"unknown negative kind {kind!r}")

    if not candidates:
        sys.exit(f"no {kind} negatives could be constructed")

    want = int(len(pairs) * ratio)
    idx = rng.permutation(len(candidates))[:want]
    neg = pd.DataFrame([candidates[i] for i in idx],
                       columns=[PEPTIDE, RECEPTOR])
    neg[FAMILY] = neg[RECEPTOR].map(fam)
    return neg, len(candidates)


def _entity_folds(col, folds, rng):
    """Assign each distinct value in `col` to a fold; return per-row folds."""
    uniq = col.unique()
    order = rng.permutation(len(uniq))
    m = {g: int(order[i] % folds) for i, g in enumerate(uniq)}
    return col.map(m)


def assign_folds(df, scheme, folds, rng):
    """Add a 'fold' column. Grouping depends on the scheme."""
    if scheme == "random_edge":
        f = rng.permutation(len(df)) % folds
        return pd.Series(f, index=df.index)

    if scheme == "cold_receptor":
        groups = df[RECEPTOR]
    elif scheme == "cold_both":
        # Partition receptors and peptides independently. A row is in test
        # fold k only if BOTH its receptor and its peptide are assigned to k;
        # it is in training only if NEITHER is. Rows where exactly one is held
        # out are discarded (fold = -1), because they would leak one entity.
        rfold = _entity_folds(df[RECEPTOR], folds, rng)
        pfold = _entity_folds(df[PEPTIDE], folds, rng)
        out = pd.Series(-1, index=df.index, dtype=int)
        for k in range(folds):
            both = (rfold == k) & (pfold == k)
            out[both] = k
        return out
    else:
        sys.exit(f"unknown scheme {scheme!r}")

    uniq = groups.unique()
    order = rng.permutation(len(uniq))
    fold_of = {g: int(order[i] % folds) for i, g in enumerate(uniq)}
    return groups.map(fold_of)


def check_leakage(df, scheme, folds, say):
    """Fail loudly if a held-out entity is also in training."""
    if scheme == "random_edge":
        return
    usable = df["fold"] >= 0
    cols = [RECEPTOR] + ([PEPTIDE] if scheme == "cold_both" else [])
    for k in range(folds):
        test_rows = df["fold"] == k
        train_rows = usable & (df["fold"] != k)
        for col in cols:
            overlap = set(df.loc[test_rows, col]) & set(df.loc[train_rows, col])
            if overlap:
                sys.exit(f"LEAKAGE in {scheme} fold {k}: {len(overlap)} "
                         f"values of {col} in both train and test")
    n_drop = int((~usable).sum())
    if n_drop:
        say(f"  {n_drop:,} rows discarded (one entity held out, not both) "
            "-- these are excluded from both train and test")
    say(f"  leakage check passed for {scheme}"
        + (f" ({', '.join(cols)})" if cols else ""))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default="data/pairs.csv")
    ap.add_argument("--out", default="data/splits")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--neg-ratio", type=float, default=1.0,
                    help="negatives per positive")
    a = ap.parse_args()

    log = []

    def say(m):
        print(m)
        log.append(m)

    pairs = load_pairs(a.pairs)
    say(f"positives: {len(pairs):,} pairs, "
        f"{pairs[RECEPTOR].nunique():,} receptors, "
        f"{pairs[PEPTIDE].nunique():,} peptides")

    n_fam = pairs[FAMILY].isna().sum()
    if n_fam:
        say(f"warning: {n_fam} pairs have no receptor_family; they cannot "
            "contribute hard negatives")

    full_grid = pairs[RECEPTOR].nunique() * pairs[PEPTIDE].nunique()
    say(f"candidate grid: {full_grid:,} possible pairs, "
        f"{full_grid - len(pairs):,} of them untested")
    say("NOTE: untested pairs are used as negatives but are not verified "
        "negatives. Some are real pairings nobody has assayed, which caps "
        "achievable AUROC.")

    os.makedirs(a.out, exist_ok=True)

    for kind in ("hard", "random"):
        for seed in a.seeds:
            rng = np.random.default_rng(seed)
            neg, pool = build_negatives(pairs, kind, a.neg_ratio, rng)
            if seed == a.seeds[0]:
                say(f"\n{kind} negatives: pool of {pool:,} candidates, "
                    f"sampling {len(neg):,} per seed "
                    f"(ratio {a.neg_ratio:g})")

            pos = pairs[[PEPTIDE, RECEPTOR, FAMILY]].copy()
            pos["label"] = 1
            neg = neg.copy()
            neg["label"] = 0
            allp = pd.concat([pos, neg], ignore_index=True)

            for scheme in ("cold_receptor", "cold_both", "random_edge"):
                rng2 = np.random.default_rng(seed + 1000)
                d = allp.copy()
                d["fold"] = assign_folds(d, scheme, a.folds, rng2)
                if seed == a.seeds[0]:
                    say(f"\n{scheme} / {kind} / seed {seed}")
                    sizes = d.groupby("fold")["label"].agg(["size", "sum"])
                    for f, row in sizes.iterrows():
                        say(f"  fold {f}: {int(row['size']):,} rows, "
                            f"{int(row['sum']):,} positive "
                            f"({row['sum']/row['size']:.3f} base rate)")
                    check_leakage(d, scheme, a.folds, say)

                fn = f"splits_{scheme}_{kind}_seed{seed}.csv"
                d.to_csv(os.path.join(a.out, fn), index=False)

    say(f"\nwrote split files to {a.out}/")
    with open(os.path.join(a.out, "splits_report.txt"), "w") as f:
        f.write("\n".join(log) + "\n")
    print("wrote", os.path.join(a.out, "splits_report.txt"))


if __name__ == "__main__":
    main()
