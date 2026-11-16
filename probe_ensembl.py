#!/usr/bin/env python3
"""
Does a second orthology source cover the proteins OMA left out?

OMA assigns no orthology group (oma_group 0, no HOG) to 19 of 248 proteins in
this dataset, concentrated in the chemokine system. This probes Ensembl
Compara for exactly those, to decide between:

  (a) add a second source and compare coverage across databases, or
  (b) exclude the 19 and document OMA coverage as a study limitation.

Scope guard: this script touches ONLY the accessions with empty OMA
fingerprints. It builds nothing and commits to nothing.

Usage:
    python probe_ensembl.py
    python probe_ensembl.py --oma data/oma_all

Note on mixing sources: Ensembl Compara spans a few hundred vertebrates while
OMA spans ~3,000 species across all domains. Fingerprints from the two are NOT
interchangeable -- any comparison has to be per-database, on a shared species
set, not merged into one matrix.

Data: Ensembl (Apache 2.0 / CC0 for data). https://rest.ensembl.org
"""

import argparse
import glob
import json
import os
import sys
import time
import urllib.error
import urllib.request

import pandas as pd

ENSEMBL = "https://rest.ensembl.org"
UA = "academic script (peptide-GPCR pairing project)"


def get(path, timeout=120):
    req = urllib.request.Request(ENSEMBL + path, headers={
        "User-Agent": UA, "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return json.loads(r.read()), None
    except urllib.error.HTTPError as e:
        return None, f"HTTP {e.code}"
    except Exception as e:
        return None, f"{type(e).__name__}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--oma", default="data/oma_all")
    ap.add_argument("--pairs", default="data/pairs_resolved.csv")
    ap.add_argument("--out", default="data")
    ap.add_argument("--pause", type=float, default=0.4)
    a = ap.parse_args()

    # which accessions have an empty OMA fingerprint
    recs = [json.load(open(f))
            for f in glob.glob(os.path.join(a.oma, "proteins", "*.json"))]
    empty = sorted({r["accession"] for r in recs
                    if r.get("status") == "ok" and not r.get("taxa")})
    if not empty:
        sys.exit(f"No empty fingerprints found under {a.oma}. Nothing to probe.")
    print(f"{len(empty)} accessions with no OMA orthology group\n")

    # map accession -> gene symbol and family from the pair table
    df = pd.read_csv(a.pairs, dtype=str)
    sym, fam = {}, {}
    for _, r in df.iterrows():
        if isinstance(r.get("receptor_uniprot"), str):
            sym.setdefault(r["receptor_uniprot"], r.get("receptor_gene"))
            fam.setdefault(r["receptor_uniprot"], r.get("receptor_family"))
        if isinstance(r.get("precursor_uniprot"), str):
            sym.setdefault(r["precursor_uniprot"], r.get("precursor_gene"))
            fam.setdefault(r["precursor_uniprot"], "(precursor)")

    rows = []
    for i, acc in enumerate(empty, 1):
        gene = sym.get(acc)
        if not isinstance(gene, str) or not gene:
            rows.append({"accession": acc, "gene": None, "family": fam.get(acc),
                         "ensembl_orthologs": None, "species": None,
                         "status": "no gene symbol"})
            print(f"  [{i:2}/{len(empty)}] {acc:<8} (no symbol)")
            continue

        data, err = get(f"/homology/symbol/human/{gene}"
                        "?type=orthologues;format=condensed"
                        ";content-type=application/json")
        if data is None:
            n, species, status = None, None, err
        else:
            homs = []
            for d in data.get("data", []):
                homs.extend(d.get("homologies", []))
            species = sorted({h.get("species") for h in homs if h.get("species")})
            n, status = len(homs), "ok"

        rows.append({"accession": acc, "gene": gene, "family": fam.get(acc),
                     "ensembl_orthologs": n,
                     "species": len(species) if species else 0,
                     "status": status})
        print(f"  [{i:2}/{len(empty)}] {acc:<8} {str(gene):<8} "
              f"{str(n):>5} orthologs  {str(len(species) if species else 0):>4} species"
              f"   {status}")
        time.sleep(a.pause)

    out = pd.DataFrame(rows)
    path = os.path.join(a.out, "ensembl_coverage_probe.csv")
    out.to_csv(path, index=False)

    got = out["species"].fillna(0).gt(0).sum()
    print(f"\ncovered by Ensembl with >0 ortholog species: {got} of {len(out)}")
    if got:
        med = out.loc[out["species"] > 0, "species"].median()
        print(f"median ortholog species among those: {int(med)}")
    print(f"\nwrote {path}")
    print("\nDecision rule: if most of these are covered, a second source is "
          "worth adding (as a separate per-database comparison, not merged "
          "into the OMA matrix). If not, exclude the 19 and document OMA "
          "coverage as a limitation.")


if __name__ == "__main__":
    main()
