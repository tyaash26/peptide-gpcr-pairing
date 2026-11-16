#!/usr/bin/env python3
"""
Score co-evolution similarity for every pair in the evaluation splits.

Produces no model and fits nothing. Each candidate (peptide, receptor) pair
gets similarity scores computed from the presence/absence fingerprints, so the
benchmark can ask whether those scores discriminate true pairs from negatives.

Metrics
-------
  bit_identity   fraction of species where the two fingerprints agree.
                 This is Foster et al.'s metric (they reported 91% for known
                 pairs against 56% for random protein-protein pairs). Shared
                 ABSENCE counts as agreement, which is why it inflates as the
                 species set grows.
  jaccard        |intersection| / |union| of present species. Ignores shared
                 absence entirely.
  pearson        correlation of the two binary vectors. Undefined when either
                 has no variance within the scope; reported as NaN.
  n_both, n_either, n_receptor, n_peptide  presence counts, so the
                 conservation confounder can be conditioned on afterwards.

Species scopes
--------------
Scores are computed several times over different column subsets, because the
choice of scope is itself a finding:
  all        every taxon in the matrix
  shared     taxa where at least one of the two genes is present
  (more scopes can be added once a taxonomy mapping exists)

Usage
-----
    python score_coevolution.py
    python score_coevolution.py --fingerprints data/fingerprints_all.csv

Outputs
-------
    data/scores_<scheme>_<negatives>_seed<N>.csv
    data/scores_report.txt
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd


def load_fingerprints(path):
    fp = pd.read_csv(path, index_col=0)
    role = fp["role"] if "role" in fp.columns else None
    mat = fp.drop(columns=[c for c in ("role",) if c in fp.columns])
    mat.columns = [str(c) for c in mat.columns]
    return mat.astype("int8"), role


def pair_scores(v_r, v_p):
    """Scores for one pair of binary vectors, over all columns given."""
    n = len(v_r)
    both = int(np.sum((v_r == 1) & (v_p == 1)))
    either = int(np.sum((v_r == 1) | (v_p == 1)))
    agree = int(np.sum(v_r == v_p))
    out = {
        "bit_identity": agree / n if n else np.nan,
        "jaccard": both / either if either else np.nan,
        "n_both": both,
        "n_either": either,
        "n_receptor": int(v_r.sum()),
        "n_peptide": int(v_p.sum()),
        "n_species": n,
    }
    if v_r.std() == 0 or v_p.std() == 0:
        out["pearson"] = np.nan
    else:
        out["pearson"] = float(np.corrcoef(v_r, v_p)[0, 1])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fingerprints", default="data/fingerprints_all.csv")
    ap.add_argument("--pairs", default="data/pairs_resolved.csv")
    ap.add_argument("--splits", default="data/splits")
    ap.add_argument("--out", default="data")
    a = ap.parse_args()

    log = []

    def say(m):
        print(m)
        log.append(m)

    mat, role = load_fingerprints(a.fingerprints)
    say(f"fingerprints: {mat.shape[0]} accessions x {mat.shape[1]} taxa "
        f"({a.fingerprints})")

    # map gtopdb ids -> accessions via the pair table
    df = pd.read_csv(a.pairs, dtype=str)
    rec_acc = dict(zip(df["receptor_gtopdb_id"], df["receptor_uniprot"]))
    pep_acc = dict(zip(df["peptide_gtopdb_id"], df["precursor_uniprot"]))

    have = set(mat.index)
    empty = {i for i in have if mat.loc[i].sum() == 0}
    say(f"accessions with an all-zero fingerprint (excluded): {len(empty)}")

    files = sorted(glob.glob(os.path.join(a.splits, "splits_*_seed*.csv")))
    if not files:
        sys.exit(f"no split files under {a.splits}")
    say(f"split files: {len(files)}")

    first = True
    for path in files:
        name = os.path.basename(path)
        m = re.match(r"splits_(.+)_(hard|random|matched)_seed(\d+)\.csv", name)
        if not m:
            say(f"  skipping unrecognised file {name}")
            continue
        scheme, negs, seed = m.groups()
        s = pd.read_csv(path, dtype=str)

        rows, dropped = [], {"no_accession": 0, "empty_fingerprint": 0}
        for _, r in s.iterrows():
            ra = rec_acc.get(r["receptor_gtopdb_id"])
            pa = pep_acc.get(r["peptide_gtopdb_id"])
            if not isinstance(ra, str) or not isinstance(pa, str) \
                    or ra not in have or pa not in have:
                dropped["no_accession"] += 1
                continue
            if ra in empty or pa in empty:
                dropped["empty_fingerprint"] += 1
                continue
            v_r = mat.loc[ra].to_numpy(dtype=float)
            v_p = mat.loc[pa].to_numpy(dtype=float)

            rec = {"receptor_gtopdb_id": r["receptor_gtopdb_id"],
                   "peptide_gtopdb_id": r["peptide_gtopdb_id"],
                   "receptor_uniprot": ra, "precursor_uniprot": pa,
                   "label": int(r["label"]), "fold": int(r["fold"])}

            for scope, cols in (("all", slice(None)),
                                ("shared", (v_r == 1) | (v_p == 1))):
                if scope == "all":
                    sr, sp = v_r, v_p
                else:
                    sr, sp = v_r[cols], v_p[cols]
                for k, val in pair_scores(sr, sp).items():
                    rec[f"{scope}__{k}"] = val
            rows.append(rec)

        out = pd.DataFrame(rows)
        dest = os.path.join(a.out, f"scores_{scheme}_{negs}_seed{seed}.csv")
        out.to_csv(dest, index=False)

        if first:
            say(f"\nexample ({name}):")
            say(f"  rows scored: {len(out):,}  dropped: {dropped}")
            if len(out):
                for metric in ("all__bit_identity", "shared__jaccard",
                               "all__pearson"):
                    g = out.groupby("label")[metric].agg(["mean", "median"])
                    say(f"  {metric}:")
                    for lab in sorted(g.index):
                        say(f"     label {lab}: mean {g.loc[lab,'mean']:.3f} "
                            f"median {g.loc[lab,'median']:.3f}")
                nan_p = out["all__pearson"].isna().sum()
                say(f"  pearson undefined for {nan_p} rows")
            first = False

    say(f"\nwrote score files to {a.out}/")
    say("")
    say("NOTE: these are scores, not results. Separation of the label means "
        "is not evidence until it survives the hard-negative split and "
        "conditioning on presence counts (n_receptor, n_peptide).")

    with open(os.path.join(a.out, "scores_report.txt"), "w") as f:
        f.write("\n".join(log) + "\n")
    print("wrote", os.path.join(a.out, "scores_report.txt"))


if __name__ == "__main__":
    main()
