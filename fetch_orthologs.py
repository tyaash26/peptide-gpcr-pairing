#!/usr/bin/env python3
"""
Fetch 1:1 orthologs from OMA for every receptor and precursor accession,
and build a presence/absence matrix.

No taxonomic filter is applied at fetch time. The full taxon set is stored so
the same data can be scored at several species scopes afterwards -- comparing
scopes is the point, and refetching is expensive.

Resumable: one JSON per accession under <cache>/proteins/. Rerun after an
interruption and it skips what is already cached.

Usage
-----
    python fetch_orthologs.py --limit 5        # try five accessions first
    python fetch_orthologs.py                  # all of them
    python fetch_orthologs.py --genomes        # also fetch the genome list

Outputs
-------
    data/oma/proteins/<ACC>.json    raw cache (gitignore this)
    data/oma/genomes.json           genome list, with --genomes
    data/fingerprints.csv           accession x taxon_id presence matrix
    data/fingerprints_report.txt

Data: OMA (Orthologous MAtrix). Altenhoff A et al. (2024)
Nucleic Acids Res 52(D1):D513-D521.
"""

import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

import pandas as pd

BASE = "https://omabrowser.org/api"
UA = "academic script (peptide-GPCR pairing project)"


def get(path, params=None, timeout=120, retries=3, pause=1.0):
    url = BASE + path
    if params:
        url += "?" + urllib.parse.urlencode(params)
    last = None
    for attempt in range(retries):
        req = urllib.request.Request(url, headers={
            "User-Agent": UA, "Accept": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read()), None
        except urllib.error.HTTPError as e:
            if e.code in (404, 403):
                return None, f"HTTP {e.code}"
            last = f"HTTP {e.code}"
        except Exception as e:
            last = f"{type(e).__name__}"
        time.sleep(pause * (attempt + 1))
    return None, last


def as_list(payload):
    """OMA returns either a bare list or a paginated dict."""
    if payload is None:
        return []
    if isinstance(payload, list):
        return payload
    if isinstance(payload, dict):
        return payload.get("results", [])
    return []


def fetch_one(acc, cache_dir, pause):
    """Resolve an accession and fetch its 1:1 orthologs. Returns a record."""
    path = os.path.join(cache_dir, f"{acc}.json")
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f), "cached"

    prot, err = get(f"/protein/{acc}/")
    if prot is None:
        rec = {"accession": acc, "status": f"unresolved ({err})",
               "omaid": None, "taxa": []}
        with open(path, "w") as f:
            json.dump(rec, f)
        return rec, "unresolved"
    time.sleep(pause)

    omaid = prot.get("omaid")
    orth, err = get(f"/protein/{omaid}/orthologs/", {"rel_type": "1:1"})
    hits = as_list(orth)
    if orth is None:
        rec = {"accession": acc, "status": f"orthologs failed ({err})",
               "omaid": omaid, "taxa": []}
    else:
        taxa, names = [], {}
        for h in hits:
            sp = h.get("species") or {}
            t = sp.get("taxon_id")
            if t is not None:
                taxa.append(int(t))
                names[str(t)] = sp.get("species")
        rec = {"accession": acc,
               "status": "ok",
               "omaid": omaid,
               "canonicalid": prot.get("canonicalid"),
               "n_orthologs": len(hits),
               "taxa": sorted(set(taxa)),
               "taxon_names": names}
    with open(path, "w") as f:
        json.dump(rec, f)
    return rec, rec["status"]


def fetch_genomes(cache_dir, pause):
    """Page through the genome list so species scope can be subset later."""
    out, page = [], 1
    while True:
        payload, err = get("/genome/", {"page": page, "per_page": 100})
        rows = as_list(payload)
        if not rows:
            if err and page == 1:
                print(f"  genome list failed: {err}")
            break
        out.extend(rows)
        print(f"  genomes page {page}: {len(rows)} (total {len(out)})")
        if len(rows) < 100:
            break
        page += 1
        time.sleep(pause)
    path = os.path.join(cache_dir, "genomes.json")
    with open(path, "w") as f:
        json.dump(out, f)
    print(f"  wrote {path} ({len(out)} genomes)")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default="data/pairs_resolved.csv")
    ap.add_argument("--cache", default="data/oma")
    ap.add_argument("--out", default="data")
    ap.add_argument("--pause", type=float, default=0.3)
    ap.add_argument("--limit", type=int, default=0,
                    help="fetch only the first N accessions (smoke test)")
    ap.add_argument("--genomes", action="store_true")
    a = ap.parse_args()

    prot_dir = os.path.join(a.cache, "proteins")
    os.makedirs(prot_dir, exist_ok=True)

    df = pd.read_csv(a.pairs, dtype=str)
    rec_acc = sorted(set(df["receptor_uniprot"].dropna()))
    pre_acc = sorted(set(df["precursor_uniprot"].dropna()))
    role = {}
    for x in rec_acc:
        role[x] = "receptor"
    for x in pre_acc:
        role.setdefault(x, "precursor")
    accs = sorted(role)
    if a.limit:
        accs = accs[:a.limit]

    print(f"{len(rec_acc)} receptor + {len(pre_acc)} precursor accessions "
          f"-> {len(role)} distinct; fetching {len(accs)}\n")

    if a.genomes:
        print("genome list:")
        fetch_genomes(a.cache, a.pause)
        print()

    records, counts = [], {}
    for i, acc in enumerate(accs, 1):
        rec, how = fetch_one(acc, prot_dir, a.pause)
        rec["role"] = role[acc]
        records.append(rec)
        counts[how.split(" ")[0]] = counts.get(how.split(" ")[0], 0) + 1
        print(f"  [{i:3}/{len(accs)}] {acc:<8} {role[acc]:<9} "
              f"{len(rec.get('taxa', [])):>4} taxa   {how}")
        if how != "cached":
            time.sleep(a.pause)

    log = [f"accessions fetched: {len(records)}",
           f"outcomes: {counts}"]

    ok = [r for r in records if r.get("status") == "ok"]
    if not ok:
        print("\nNothing resolved. Stop and inspect before building a matrix.")
        return

    all_taxa = sorted({t for r in ok for t in r["taxa"]})
    log.append(f"distinct taxa seen across all accessions: {len(all_taxa):,}")
    n = pd.Series([len(r["taxa"]) for r in ok])
    log.append(f"orthologs per accession: median {int(n.median())}, "
               f"min {int(n.min())}, max {int(n.max())}")

    # presence/absence matrix: rows = accessions, columns = taxon ids
    mat = pd.DataFrame(0, index=[r["accession"] for r in ok],
                       columns=all_taxa, dtype="int8")
    for r in ok:
        mat.loc[r["accession"], r["taxa"]] = 1
    mat.insert(0, "role", [r["role"] for r in ok])

    path = os.path.join(a.out, "fingerprints.csv")
    mat.to_csv(path)
    log.append(f"matrix: {mat.shape[0]} accessions x {len(all_taxa)} taxa")
    log.append("")
    log.append("NOTE: no taxonomic filter was applied. Score at several "
               "species scopes from this one matrix; adding taxa where both "
               "genes are absent inflates bit-identity similarity.")

    print("\n" + "\n".join(log))
    print(f"\nwrote {path}")
    with open(os.path.join(a.out, "fingerprints_report.txt"), "w") as f:
        f.write("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
