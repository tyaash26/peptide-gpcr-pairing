#!/usr/bin/env python3
"""
Re-run the benchmark split by precursor type.

Tests one hypothesis: that the conservation confounder (receptor conservation
alone scoring AUROC 0.629) is driven by moonlighting proteins -- housekeeping
and structural genes whose fragments act at GPCRs, which are conserved for
reasons unrelated to their receptor.

Network analysis says those are only 35 of 386 pairs (9%), so the hypothesis is
probably wrong. This measures it rather than assuming.

Strata:
    all          every scored pair
    classical    precursor is a classical prohormone
    moonlighting precursor is in the hand-curated moonlighting set

Usage:
    python benchmark_strata.py
    python benchmark_strata.py --scheme cold_receptor --negatives hard

Outputs:
    data/benchmark_strata.csv
    data/benchmark_strata_report.txt
"""

import argparse
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

# same set as analyze_network.py -- hand-curated from gene names, not an
# annotation source. A write-up using this split needs a defensible criterion.
MOONLIGHTING = {"HBB", "PSAP", "RPS19", "YARS1", "ANXA1", "PRNP", "MT-RNR2",
                "WNT3A", "WNT4", "WNT5A", "RSPO1", "RSPO2", "RSPO3", "RSPO4",
                "C3", "C5", "SAA1", "CAMP", "HEBP1"}

METRICS = ["all__bit_identity", "shared__jaccard", "all__pearson",
           "shared__pearson", "all__n_receptor", "all__n_peptide"]


def load_benchmark_module(path="benchmark.py"):
    if not os.path.exists(path):
        sys.exit(f"{path} not found; run from the project root")
    spec = importlib.util.spec_from_file_location("bench", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", default="data")
    ap.add_argument("--pairs", default="data/pairs_resolved.csv")
    ap.add_argument("--scheme", default="cold_receptor")
    ap.add_argument("--negatives", default="hard")
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--out", default="data")
    a = ap.parse_args()

    b = load_benchmark_module()

    log = []

    def say(m=""):
        print(m)
        log.append(m)

    # accession -> precursor gene, so score files can be stratified
    pairs = pd.read_csv(a.pairs, dtype=str)
    gene_of = dict(zip(pairs["precursor_uniprot"], pairs["precursor_gene"]))
    moon_acc = {acc for acc, g in gene_of.items() if g in MOONLIGHTING}
    say(f"moonlighting precursor accessions: {len(moon_acc)} "
        f"of {len(gene_of)}")

    rows = []
    for seed in a.seeds:
        path = os.path.join(
            a.scores, f"scores_{a.scheme}_{a.negatives}_seed{seed}.csv")
        if not os.path.exists(path):
            say(f"missing {path}, skipping")
            continue
        d = pd.read_csv(path)
        d["moon"] = d["precursor_uniprot"].isin(moon_acc)

        for stratum, sub in (("all", d),
                             ("classical", d[~d.moon]),
                             ("moonlighting", d[d.moon])):
            if sub.empty or sub["label"].nunique() < 2:
                if seed == a.seeds[0]:
                    say(f"  {stratum}: too few rows or one class only "
                        f"(n={len(sub)}) -- cannot evaluate")
                continue
            y = sub["label"].to_numpy()
            for metric in METRICS:
                if metric not in sub.columns:
                    continue
                s = sub[metric].to_numpy(dtype=float)
                ro = b.auroc(y, s)
                lo, hi = b.boot_ci(y, s, b.auroc, n=a.boot, seed=seed)
                rows.append({"stratum": stratum, "seed": seed,
                             "metric": metric, "n": len(sub),
                             "n_pos": int((y == 1).sum()),
                             "auroc": ro, "lo": lo, "hi": hi,
                             "auprc": b.auprc(y, s)})

    if not rows:
        sys.exit("nothing evaluated")

    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(a.out, "benchmark_strata.csv"), index=False)

    say()
    say(f"{a.scheme} / {a.negatives} negatives, mean over "
        f"{res.seed.nunique()} seed(s)")
    for stratum, g in res.groupby("stratum", sort=False):
        n = int(g["n"].mean())
        npos = int(g["n_pos"].mean())
        say()
        say(f"{stratum}  (n={n}, positives={npos})")
        for metric, gm in g.groupby("metric"):
            mu = gm.auroc.mean()
            lo, hi = gm.lo.mean(), gm.hi.mean()
            tag = "  [confounder]" if metric.endswith(("n_receptor",
                                                       "n_peptide")) else ""
            say(f"   {metric:<22} AUROC {mu:.3f} "
                f"[{lo:.3f}–{hi:.3f}]  AUPRC {gm.auprc.mean():.3f}{tag}")

    say()
    say("Reading this. If the confounder's AUROC stays high in the "
        "'classical' stratum, moonlighting proteins are not driving it and "
        "the conservation effect is a property of the whole gene set. The "
        "'moonlighting' stratum will be small, so treat its interval as "
        "indicative only.")

    with open(os.path.join(a.out, "benchmark_strata_report.txt"), "w") as f:
        f.write("\n".join(log) + "\n")
    print(f"\nwrote {a.out}/benchmark_strata.csv and "
          f"{a.out}/benchmark_strata_report.txt")


if __name__ == "__main__":
    main()
