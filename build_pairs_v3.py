#!/usr/bin/env python3
"""
Milestone 1: human peptide-GPCR pair table from GtoPdb bulk files.  (v3)

A ligand counts as a peptide if its Ligand ID appears in peptides.csv.
('Ligand Type' is empty in gpcr_interactions.csv, so it cannot be used.)

v3 changes:
  * peptides.csv 'Species' is COMMA-delimited while 'UniProt id' is
    PIPE-delimited, and the accession order is not reliably species-aligned
    (angiotensin I-III: P01015|P11859|P01019; angiotensin IV: P11859|P01015|
    P01019 - same precursor, reordered).  So no human accession is asserted.
    The raw string is kept as precursor_uniprot_raw and resolution is deferred
    to the UniProt ID-mapping service, keyed on peptide_gene.
  * gpcr_interactions.csv 'Ligand Species' is PIPE-delimited and IS aligned
    with 'Ligand Gene Symbol', so peptide_gene is taken from there.
  * rows whose Ligand Species is blank are kept and flagged species_assumed.

Inputs (in --data, default data/gtopdb/):
    gpcr_interactions.csv, peptides.csv, targets_and_families.csv
Outputs (in --out, default data/):
    pairs.csv, pairs_report.txt

Run:  python build_pairs_v3.py
      python build_pairs_v3.py --no-endogenous-only

Data: IUPHAR/BPS Guide to PHARMACOLOGY (GtoPdb), ODbL / CC-BY-SA 4.0.
Reference: Harding SD et al. (2026) Nucleic Acids Res 54(D1):D1446-D1456.
"""

import argparse
import os
import re
import sys

import pandas as pd


def read_gtopdb(path, first_col):
    with open(path, encoding="utf-8", errors="replace") as f:
        for i in range(5):
            if f.readline().lstrip().startswith('"' + first_col + '"'):
                return pd.read_csv(path, skiprows=i, dtype=str)
    sys.exit(f"{path}: no header row starting with {first_col!r}")


def split_any(s):
    """GtoPdb uses '|' in some columns and ', ' in others."""
    if not isinstance(s, str):
        return []
    return [x.strip() for x in re.split(r"[|,]", s) if x.strip()]


def human_aligned(species, values):
    """Element of `values` aligned to 'Human' in `species`.

    Returns (value, how) where how is 'single', 'aligned' or None.
    No fallback guessing: if the lists do not line up, we return nothing
    rather than assert a wrong identifier.
    """
    sp, vals = split_any(species), split_any(values)
    if len(vals) == 1 and len(sp) <= 1:
        return vals[0], "single"
    if len(vals) == 1:
        return vals[0], "single"
    if len(sp) == len(vals):
        for s, v in zip(sp, vals):
            if s.lower() == "human":
                return v, "aligned"
        return None, None
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/gtopdb")
    ap.add_argument("--out", default="data")
    ap.add_argument("--no-endogenous-only", action="store_true")
    a = ap.parse_args()

    log = []

    def say(m):
        print(m)
        log.append(m)

    inter = read_gtopdb(os.path.join(a.data, "gpcr_interactions.csv"), "Target")
    pep = read_gtopdb(os.path.join(a.data, "peptides.csv"), "Ligand id")
    tgt = read_gtopdb(os.path.join(a.data, "targets_and_families.csv"), "Type")
    say(f"gpcr_interactions: {len(inter):,} rows, {inter.shape[1]} cols")
    say(f"peptides:          {len(pep):,} rows")
    say(f"targets_families:  {len(tgt):,} rows")

    pep = pep.drop_duplicates(subset="Ligand id")
    pep_ids = set(pep["Ligand id"].dropna())
    say(f"\ndistinct peptide ligand IDs:              {len(pep_ids):,}")

    df = inter[inter["Ligand ID"].isin(pep_ids)]
    say(f"GPCR interactions with a peptide ligand:  {len(df):,}")

    df = df[df["Target Species"].fillna("").str.strip().eq("Human")]
    say(f"  + human receptor:                       {len(df):,}")

    ls = df["Ligand Species"].fillna("")
    keep = ls.str.lower().str.contains("human") | ls.str.strip().eq("")
    say(f"  + human peptide (blanks kept, flagged): {int(keep.sum()):,}"
        f"   [{int(ls.str.strip().eq('').sum()):,} blank]")
    df = df[keep]

    if not a.no_endogenous_only:
        df = df[df["Endogenous"].fillna("").str.strip().str.lower().eq("true")]
        say(f"  + endogenous == true:                   {len(df):,}")

    if df.empty:
        sys.exit("No rows survived. See counts above.")

    say("\nType:   " + str(dict(df["Type"].fillna("(blank)").value_counts())))
    say("Action: " + str(dict(df["Action"].fillna("(blank)").value_counts())))

    df = df.copy()
    gene = [human_aligned(s, g) for s, g in
            zip(df["Ligand Species"], df["Ligand Gene Symbol"])]
    df["peptide_gene"] = [g for g, _ in gene]
    df["peptide_gene_basis"] = [h for _, h in gene]
    df["species_assumed"] = df["Ligand Species"].fillna("").str.strip().eq("")

    p = pep[["Ligand id", "Species", "Single letter amino acid sequence",
             "UniProt id", "Ensembl id",
             "Post-translational modification"]].rename(
        columns={"Ligand id": "Ligand ID",
                 "Species": "pep_species",
                 "Single letter amino acid sequence": "peptide_sequence",
                 "UniProt id": "precursor_uniprot_raw",
                 "Ensembl id": "precursor_ensembl_raw",
                 "Post-translational modification": "peptide_ptm"})
    df = df.merge(p, on="Ligand ID", how="left")

    t = (tgt[["Target id", "Family name", "HGNC symbol"]]
         .rename(columns={"Target id": "Target ID",
                          "Family name": "receptor_family",
                          "HGNC symbol": "receptor_hgnc"})
         .drop_duplicates(subset="Target ID"))
    df = df.merge(t, on="Target ID", how="left")

    out = pd.DataFrame({
        "receptor_name": df["Target"],
        "receptor_gtopdb_id": df["Target ID"],
        "receptor_gene": df["Target Gene Symbol"].fillna(df["receptor_hgnc"]),
        "receptor_uniprot": df["Target UniProt ID"],
        "receptor_ensembl": df["Target Ensembl Gene ID"],
        "receptor_family": df["receptor_family"],
        "peptide_name": df["Ligand"],
        "peptide_gtopdb_id": df["Ligand ID"],
        "peptide_gene": df["peptide_gene"],
        "peptide_gene_basis": df["peptide_gene_basis"],
        "peptide_sequence": df["peptide_sequence"],
        "precursor_uniprot_raw": df["precursor_uniprot_raw"],
        "precursor_ensembl_raw": df["precursor_ensembl_raw"],
        "peptide_ptm": df["peptide_ptm"],
        "action": df["Action"],
        "type": df["Type"],
        "endogenous": df["Endogenous"],
        "primary_target": df["Primary Target"],
        "affinity_median": df["Affinity Median"],
        "affinity_units": df["Affinity Units"],
        "pubmed_id": df["PubMed ID"],
        "species_assumed": df["species_assumed"],
    })

    before = len(out)
    out = (out.sort_values("affinity_median", na_position="last")
              .drop_duplicates(subset=["peptide_gtopdb_id", "receptor_gtopdb_id"])
              .sort_values(["receptor_gene", "peptide_name"]))
    say(f"\n{before:,} rows -> {len(out):,} unique peptide-receptor pairs")
    say(f"  distinct receptors:  {out['receptor_gtopdb_id'].nunique():,}")
    say(f"  distinct peptides:   {out['peptide_gtopdb_id'].nunique():,}")
    say(f"  no sequence:         {out['peptide_sequence'].isna().sum():,}")
    say(f"  no precursor gene:   {out['peptide_gene'].isna().sum():,}"
        "   (blocks milestone 2 for those)")
    say(f"  species_assumed:     {int(out['species_assumed'].sum()):,}")
    say("  peptide_gene basis:  "
        + str(dict(out["peptide_gene_basis"].fillna("(none)").value_counts())))
    say("\nNOTE: precursor_uniprot_raw is the unparsed GtoPdb string. The "
        "accession order is not species-aligned, so no human accession is "
        "asserted here; resolve from peptide_gene via UniProt ID mapping.")

    os.makedirs(a.out, exist_ok=True)
    out.to_csv(os.path.join(a.out, "pairs.csv"), index=False)
    with open(os.path.join(a.out, "pairs_report.txt"), "w") as f:
        f.write("\n".join(log) + "\n")
    say(f"\nwrote {a.out}/pairs.csv and {a.out}/pairs_report.txt")

    print("\n--- angiotensin check (should show no accession) ---")
    chk = out[out["peptide_name"].str.contains("angiotensin", case=False, na=False)]
    print(chk[["receptor_gene", "peptide_name", "peptide_gene",
               "precursor_uniprot_raw"]].head(8).to_string(index=False))


if __name__ == "__main__":
    main()
