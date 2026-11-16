#!/usr/bin/env python3
"""
Resolve peptide precursor gene symbols to reviewed human UniProt accessions.

build_pairs_v3.py deliberately asserts no precursor accession, because GtoPdb
lists them in an order that is not species-aligned. This script resolves them
properly from the gene symbol via the UniProt REST API, then CROSS-CHECKS each
answer against the raw GtoPdb string that was preserved in pairs.csv: if
UniProt returns P01019 for AGT and P01019 appears among the raw pipe-delimited
accessions, the two sources agree independently.

Outputs (in --out, default data/):
    precursors.csv        gene, accession, agrees_with_gtopdb, status
    precursors_report.txt

Usage:
    python resolve_precursors.py
    python resolve_precursors.py --pause 0.3

Data: UniProtKB (CC-BY 4.0) and GtoPdb (ODbL / CC-BY-SA 4.0).
"""

import argparse
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

import pandas as pd

UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
UA = "academic script (peptide-GPCR pairing project)"


def query_uniprot(gene, timeout=60):
    """Reviewed human entries whose primary gene name is exactly `gene`."""
    params = {
        "query": f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true",
        "fields": "accession,gene_primary,protein_name,length",
        "format": "json",
        "size": "10",
    }
    url = UNIPROT + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={
        "User-Agent": UA, "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            data = json.loads(r.read())
    except urllib.error.HTTPError as e:
        return None, f"HTTP {e.code}"
    except Exception as e:
        return None, f"{type(e).__name__}"
    results = data.get("results", [])
    if not results:
        return None, "no reviewed human entry"
    hits = []
    for r in results:
        acc = r.get("primaryAccession")
        primary = ""
        for g in r.get("genes", []) or []:
            primary = (g.get("geneName") or {}).get("value", "")
            break
        if primary.upper() == gene.upper():
            hits.append(acc)
    if not hits:
        return None, f"{len(results)} hits, none with primary name {gene}"
    if len(hits) > 1:
        return hits[0], f"ambiguous ({len(hits)} reviewed entries)"
    return hits[0], "ok"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default="data/pairs.csv")
    ap.add_argument("--out", default="data")
    ap.add_argument("--pause", type=float, default=0.3)
    a = ap.parse_args()

    df = pd.read_csv(a.pairs, dtype=str)
    for c in ("peptide_gene", "precursor_uniprot_raw"):
        if c not in df.columns:
            sys.exit(f"{a.pairs}: missing column {c!r}")

    # one row per precursor gene, keeping every raw string seen for it
    raw_by_gene = {}
    for g, raw in zip(df["peptide_gene"], df["precursor_uniprot_raw"]):
        if not isinstance(g, str) or not g.strip():
            continue
        s = raw_by_gene.setdefault(g.strip(), set())
        if isinstance(raw, str):
            for v in re.split(r"[|,]", raw):
                v = v.strip()
                if v:
                    s.add(v)

    genes = sorted(raw_by_gene)
    print(f"{len(genes)} distinct precursor genes to resolve\n")

    rows, log = [], [f"{len(genes)} distinct precursor genes"]
    agree = disagree = absent = failed = 0

    for i, gene in enumerate(genes, 1):
        acc, status = query_uniprot(gene)
        raws = raw_by_gene[gene]
        if acc is None:
            check = "no_accession"
            failed += 1
        elif not raws:
            check = "no_raw_to_check"
            absent += 1
        elif acc in raws:
            check = "agrees"
            agree += 1
        else:
            check = "DISAGREES"
            disagree += 1
        rows.append({"peptide_gene": gene, "precursor_uniprot": acc,
                     "status": status, "gtopdb_check": check,
                     "gtopdb_raw_values": "|".join(sorted(raws))})
        flag = "  <-- check" if check == "DISAGREES" else ""
        print(f"  [{i:3}/{len(genes)}] {gene:<12} {str(acc):<10} "
              f"{check:<15} {status}{flag}")
        time.sleep(a.pause)

    out = pd.DataFrame(rows)
    os.makedirs(a.out, exist_ok=True)
    out.to_csv(os.path.join(a.out, "precursors.csv"), index=False)

    log += [
        f"agrees with a GtoPdb raw accession: {agree}",
        f"disagrees with every GtoPdb raw accession: {disagree}",
        f"no raw string available to check: {absent}",
        f"unresolved by UniProt: {failed}",
        "",
        "A 'DISAGREES' row means UniProt's reviewed human accession for that "
        "gene is absent from the accessions GtoPdb listed. Inspect these "
        "individually; do not auto-correct either source.",
    ]
    print("\n" + "\n".join(log[1:]))
    with open(os.path.join(a.out, "precursors_report.txt"), "w") as f:
        f.write("\n".join(log) + "\n")
    print(f"\nwrote {a.out}/precursors.csv and {a.out}/precursors_report.txt")


if __name__ == "__main__":
    main()
