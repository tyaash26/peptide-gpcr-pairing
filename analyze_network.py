#!/usr/bin/env python3
"""
Describe the structure of the known peptide-GPCR network.

No orthology, no model, no network calls: this reads pairs_resolved.csv and
reports how the bipartite graph is shaped. It answers questions the benchmark
cannot -- how promiscuous peptides and receptors are, which families are
one-to-one, and whether the moonlighting proteins behave differently.

Usage:
    python analyze_network.py

Outputs (in --out, default data/):
    network_degrees.csv     per-entity degree, with family
    network_report.txt
"""

import argparse
import os

import pandas as pd

# Entries that are not classical prohormones: full proteins whose fragments
# act at GPCRs, or Frizzled/LGR ligands. Flagged because they are conserved
# for reasons unrelated to their receptor.
MOONLIGHTING = {"HBB", "PSAP", "RPS19", "YARS1", "ANXA1", "PRNP", "MT-RNR2",
                "WNT3A", "WNT4", "WNT5A", "RSPO1", "RSPO2", "RSPO3", "RSPO4",
                "C3", "C5", "SAA1", "CAMP", "HEBP1"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default="data/pairs_resolved.csv")
    ap.add_argument("--out", default="data")
    a = ap.parse_args()

    log = []

    def say(m=""):
        print(m)
        log.append(m)

    df = pd.read_csv(a.pairs, dtype=str)
    say(f"{len(df):,} pairs · "
        f"{df['receptor_gtopdb_id'].nunique():,} receptors · "
        f"{df['peptide_gtopdb_id'].nunique():,} peptides · "
        f"{df['precursor_uniprot'].nunique():,} precursor genes")

    # --- degree distributions -------------------------------------------
    rec_deg = df.groupby("receptor_gtopdb_id").agg(
        name=("receptor_gene", "first"),
        family=("receptor_family", "first"),
        degree=("peptide_gtopdb_id", "nunique")).reset_index()
    pep_deg = df.groupby("peptide_gtopdb_id").agg(
        name=("peptide_name", "first"),
        precursor=("precursor_gene", "first"),
        degree=("receptor_gtopdb_id", "nunique")).reset_index()

    say()
    say("Receptors per peptide")
    say(f"  median {pep_deg.degree.median():.0f}, "
        f"mean {pep_deg.degree.mean():.2f}, max {pep_deg.degree.max()}")
    vc = pep_deg.degree.value_counts().sort_index()
    for k, v in vc.items():
        say(f"    {k} receptor(s): {v:>3} peptides "
            f"({v/len(pep_deg)*100:.0f}%)")

    say()
    say("Peptides per receptor")
    say(f"  median {rec_deg.degree.median():.0f}, "
        f"mean {rec_deg.degree.mean():.2f}, max {rec_deg.degree.max()}")
    vc = rec_deg.degree.value_counts().sort_index()
    for k, v in vc.items():
        say(f"    {k} peptide(s): {v:>3} receptors "
            f"({v/len(rec_deg)*100:.0f}%)")

    say()
    say("Most promiscuous peptides")
    for _, r in pep_deg.nlargest(8, "degree").iterrows():
        say(f"  {r['name'][:34]:<34} {r.degree:>2} receptors "
            f"(from {r.precursor})")

    say()
    say("Most promiscuous receptors")
    for _, r in rec_deg.nlargest(8, "degree").iterrows():
        say(f"  {str(r['name']):<10} {r.degree:>2} peptides "
            f"({str(r.family)[:40]})")

    # --- by family -------------------------------------------------------
    say()
    say("By receptor family (families with >= 3 receptors)")
    fam = df.groupby("receptor_family").agg(
        receptors=("receptor_gtopdb_id", "nunique"),
        peptides=("peptide_gtopdb_id", "nunique"),
        pairs=("receptor_gtopdb_id", "size")).reset_index()
    fam["pairs_per_receptor"] = (fam.pairs / fam.receptors).round(2)
    fam = fam[fam.receptors >= 3].sort_values("pairs_per_receptor",
                                              ascending=False)
    for _, r in fam.iterrows():
        say(f"  {str(r['receptor_family'])[:38]:<38} "
            f"{r.receptors:>3} rec {r.peptides:>3} pep "
            f"{r.pairs:>3} pairs  {r.pairs_per_receptor:>5} per receptor")

    # --- precursors yielding several peptides ----------------------------
    say()
    say("Precursors yielding several distinct peptides")
    multi = (df.groupby("precursor_gene")["peptide_gtopdb_id"].nunique()
               .sort_values(ascending=False))
    multi = multi[multi > 1]
    say(f"  {len(multi)} of {df['precursor_gene'].nunique()} precursors "
        "give more than one peptide")
    for g, n in multi.head(8).items():
        say(f"    {g:<10} {n} peptides")

    # --- moonlighting vs classical ---------------------------------------
    say()
    say("Classical prohormones against moonlighting proteins")
    df["is_moonlighting"] = df["precursor_gene"].isin(MOONLIGHTING)
    for flag, label in ((False, "classical prohormone"),
                        (True, "moonlighting protein")):
        sub = df[df.is_moonlighting == flag]
        if sub.empty:
            continue
        say(f"  {label:<22} {len(sub):>3} pairs, "
            f"{sub['precursor_gene'].nunique():>3} precursors, "
            f"{sub['receptor_gtopdb_id'].nunique():>3} receptors")
    say("  (moonlighting set is hand-curated; see MOONLIGHTING in this script)")

    # --- write -----------------------------------------------------------
    os.makedirs(a.out, exist_ok=True)
    rec_deg["side"] = "receptor"
    pep_deg["side"] = "peptide"
    out = pd.concat([
        rec_deg.rename(columns={"receptor_gtopdb_id": "gtopdb_id"}),
        pep_deg.rename(columns={"peptide_gtopdb_id": "gtopdb_id"}),
    ], ignore_index=True)
    out.to_csv(os.path.join(a.out, "network_degrees.csv"), index=False)
    with open(os.path.join(a.out, "network_report.txt"), "w") as f:
        f.write("\n".join(log) + "\n")
    print(f"\nwrote {a.out}/network_degrees.csv and {a.out}/network_report.txt")


if __name__ == "__main__":
    main()
