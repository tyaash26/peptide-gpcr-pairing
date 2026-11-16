#!/usr/bin/env python3
"""
Final data-assembly step: key the pair table on UniProt accessions.

Why this exists
---------------
GtoPdb gene symbols cannot be trusted for identity in this file:

  * 'Ligand Type' is empty throughout gpcr_interactions.csv
  * peptides.csv 'Species' is comma-delimited while 'UniProt id' is
    pipe-delimited, so the two never align
  * 'Ligand Gene Symbol' is pipe-delimited and the SAME LENGTH as
    'Ligand Species', but not in the same order -- neurotensin (ligand 1579)
    has species Bovine|Human|Mouse|Rat against symbols NTS|Nts|Nts|NTS
  * ligand 1866 (prokineticin-1) carries prokineticin-2's accession

Symbol matching is therefore retired. Accessions resolved from UniProt and
corroborated against the preserved raw GtoPdb string (see
resolve_precursors.py) become the join key, which collapses NTS/Nts
automatically and is also what the OMA API resolves by.

Inputs (in --data, default data/):
    pairs.csv, precursors.csv, overrides.csv (optional)
Outputs:
    pairs_resolved.csv, pairs_resolved_report.txt

Usage:
    python resolve_pairs.py
"""

import argparse
import os
import sys

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--overrides", default="overrides.csv")
    a = ap.parse_args()

    log = []

    def say(m):
        print(m)
        log.append(m)

    pairs = pd.read_csv(os.path.join(a.data, "pairs.csv"), dtype=str)
    prec = pd.read_csv(os.path.join(a.data, "precursors.csv"), dtype=str)
    say(f"pairs.csv:      {len(pairs):,} rows")
    say(f"precursors.csv: {len(prec):,} genes")

    # --- attach the resolved precursor accession ------------------------
    m = prec[["peptide_gene", "precursor_uniprot", "gtopdb_check", "status"]]
    df = pairs.merge(m, on="peptide_gene", how="left")

    # --- apply documented overrides -------------------------------------
    n_over = 0
    if os.path.exists(a.overrides):
        ov = pd.read_csv(a.overrides, dtype=str)
        say(f"\noverrides.csv: {len(ov)} entries")
        for _, r in ov.iterrows():
            mask = df["peptide_gtopdb_id"].eq(r["peptide_gtopdb_id"])
            if not mask.any():
                say(f"  warning: override for ligand {r['peptide_gtopdb_id']} "
                    "matched no rows")
                continue
            before = df.loc[mask, r["field"]].dropna().unique()
            df.loc[mask, r["field"]] = r["value"]
            if r["field"] == "precursor_uniprot":
                # the GtoPdb cross-check referred to the superseded value
                df.loc[mask, "gtopdb_check"] = "overridden"
            n_over += int(mask.sum())
            say(f"  ligand {r['peptide_gtopdb_id']}: {r['field']} "
                f"{list(before)} -> {r['value']} ({int(mask.sum())} rows)")
    else:
        say(f"\nno {a.overrides} found; proceeding without overrides")

    # --- symbol collapse, verified by accession -------------------------
    no_acc = df["precursor_uniprot"].isna().sum()
    say(f"\npairs with no precursor accession: {no_acc:,}")

    sym_per_acc = (df.dropna(subset=["precursor_uniprot"])
                     .groupby("precursor_uniprot")["peptide_gene"]
                     .nunique())
    collapsed = sym_per_acc[sym_per_acc > 1]
    if len(collapsed):
        say(f"accessions carrying more than one symbol "
            f"(collapsed by this step): {len(collapsed)}")
        for acc in collapsed.index:
            syms = sorted(set(df.loc[df["precursor_uniprot"].eq(acc),
                                     "peptide_gene"].dropna()))
            say(f"  {acc}: {syms}")
    else:
        say("no accession carries more than one symbol")

    acc_per_sym = (df.dropna(subset=["precursor_uniprot"])
                     .groupby("peptide_gene")["precursor_uniprot"].nunique())
    split = acc_per_sym[acc_per_sym > 1]
    if len(split):
        say(f"WARNING: {len(split)} symbols map to several accessions: "
            f"{list(split.index)}")

    # canonical display symbol: the uppercase form seen for that accession
    def canon(group):
        syms = [s for s in group.dropna().unique()]
        upper = [s for s in syms if s.isupper()]
        return (upper or syms or [None])[0]

    canon_map = (df.dropna(subset=["precursor_uniprot"])
                   .groupby("precursor_uniprot")["peptide_gene"]
                   .apply(canon).to_dict())
    df["precursor_gene"] = df["precursor_uniprot"].map(canon_map)

    # --- deduplicate on (peptide, receptor) -----------------------------
    # NOT on (precursor, receptor): one precursor yields several distinct
    # peptides with different pharmacology -- GCG gives glucagon and GLP-1,
    # which act at different receptors; PDYN gives several dynorphins; GAST,
    # CXCL12 and TAC1 likewise. Keying on the precursor merged them and lost
    # 99 real pairs. The precursor accession is a FEATURE key, carried
    # alongside; the pair is peptide-to-receptor.
    key = ["peptide_gtopdb_id", "receptor_uniprot"]
    have = df.dropna(subset=key)
    before = len(df)
    dedup = (have.sort_values("affinity_median", na_position="last")
                 .drop_duplicates(subset=key))
    lost = before - len(df.dropna(subset=key))
    say(f"\n{before:,} rows -> {len(dedup):,} unique "
        "(peptide, receptor) pairs")
    if lost:
        say(f"  {lost:,} rows dropped for want of an accession on one side")

    say(f"  distinct receptors:  {dedup['receptor_uniprot'].nunique():,}")
    say(f"  distinct peptides:   {dedup['peptide_gtopdb_id'].nunique():,}")
    say(f"  distinct precursors: {dedup['precursor_uniprot'].nunique():,}"
        "  (feature key; several peptides may share one)")
    say(f"  rows carrying an override: {n_over:,}")

    disagree = dedup["gtopdb_check"].eq("DISAGREES").sum()
    say(f"  pairs whose precursor still disagrees with GtoPdb: {disagree:,}"
        " (should be 0 once overridden)")

    cols = ["receptor_name", "receptor_gtopdb_id", "receptor_gene",
            "receptor_uniprot", "receptor_ensembl", "receptor_family",
            "peptide_name", "peptide_gtopdb_id", "precursor_gene",
            "precursor_uniprot", "peptide_sequence", "peptide_ptm",
            "action", "type", "endogenous", "primary_target",
            "affinity_median", "affinity_units", "pubmed_id"]
    out = dedup[[c for c in cols if c in dedup.columns]].sort_values(
        ["receptor_gene", "peptide_name"])

    path = os.path.join(a.data, "pairs_resolved.csv")
    out.to_csv(path, index=False)
    say(f"\nwrote {path}")

    with open(os.path.join(a.data, "pairs_resolved_report.txt"), "w") as f:
        f.write("\n".join(log) + "\n")
    print("wrote", os.path.join(a.data, "pairs_resolved_report.txt"))

    print("\n--- neurotensin check (NTS/Nts should now be one gene) ---")
    chk = out[out["peptide_name"].str.contains("neurotensin", case=False,
                                               na=False)]
    print(chk[["receptor_gene", "peptide_name", "precursor_gene",
               "precursor_uniprot"]].to_string(index=False))


if __name__ == "__main__":
    main()
