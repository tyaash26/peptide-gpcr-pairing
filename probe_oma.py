#!/usr/bin/env python3
"""
Probe the OMA REST API before writing the ortholog fetcher.

Resolves nothing and assumes nothing: it tries several documented identifier
forms for one known human gene and reports which return usable JSON, plus the
shape of the ortholog response and the species coverage.

Test gene: BDKRB2 / UniProt P30411 (a peptide GPCR with a known partner).

Usage:
    python probe_oma.py
    python probe_oma.py --gene AGTR1 --uniprot P30556

Data: OMA (Orthologous MAtrix). Altenhoff A et al. (2024)
Nucleic Acids Res 52(D1):D513-D521.
"""

import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

BASE = "https://omabrowser.org/api"
UA = "academic script (peptide-GPCR pairing project)"


def get(path, params=None, timeout=60):
    url = BASE + path
    if params:
        url += "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={
        "User-Agent": UA, "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            code = r.status
            raw = r.read()
    except urllib.error.HTTPError as e:
        return e.code, None, url
    except Exception as e:
        return None, f"{type(e).__name__}: {e}", url
    try:
        return code, json.loads(raw), url
    except json.JSONDecodeError:
        return code, f"non-JSON ({len(raw)} bytes)", url


def show(label, code, data, url):
    if isinstance(data, str) or data is None:
        print(f"  [{str(code):>4}] {label:<34} {data or ''}")
        return False
    if isinstance(data, dict):
        keys = list(data)[:8]
        n = data.get("count", "")
        print(f"  [{code:>4}] {label:<34} dict keys={keys}"
              + (f" count={n}" if n != "" else ""))
    elif isinstance(data, list):
        print(f"  [{code:>4}] {label:<34} list of {len(data)}")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", default="BDKRB2")
    ap.add_argument("--uniprot", default="P30411")
    ap.add_argument("--pause", type=float, default=0.5)
    a = ap.parse_args()

    print(f"probing {BASE}  (gene {a.gene}, uniprot {a.uniprot})\n")

    print("1. identifier resolution")
    attempts = [
        (f"/protein/{a.uniprot}/",            "protein by UniProt"),
        (f"/protein/{a.gene}_HUMAN/",         "protein by GENE_HUMAN"),
        (f"/protein/HUMAN{a.gene}/",          "protein by HUMAN+symbol"),
    ]
    resolved = None
    for path, label in attempts:
        code, data, url = get(path)
        ok = show(label, code, data, url)
        if ok and isinstance(data, dict) and resolved is None:
            resolved = data
        time.sleep(a.pause)

    for path, label, params in [
        ("/xref/", "xref search", {"search": a.gene}),
        ("/search/", "generic search", {"query": a.gene}),
    ]:
        code, data, url = get(path, params)
        show(label, code, data, url)
        time.sleep(a.pause)

    if resolved is None:
        print("\nNo identifier form resolved. Stop here and report the codes "
              "above rather than guessing further.")
        return

    omaid = resolved.get("omaid") or resolved.get("entry_nr")
    print(f"\nresolved entry: omaid={resolved.get('omaid')} "
          f"entry_nr={resolved.get('entry_nr')} "
          f"canonicalid={resolved.get('canonicalid')}")
    print(f"  species: {(resolved.get('species') or {}).get('species')}")

    print("\n2. ortholog endpoints")
    for path, label, params in [
        (f"/protein/{omaid}/orthologs/", "pairwise orthologs", None),
        (f"/protein/{omaid}/orthologs/", "orthologs rel_type 1:1",
         {"rel_type": "1:1"}),
        (f"/protein/{omaid}/hog_derived_orthologs/", "hog-derived orthologs",
         {"per_page": 5}),
    ]:
        code, data, url = get(path, params)
        show(label, code, data, url)
        if isinstance(data, dict) and "results" in data:
            res = data["results"]
            if res:
                print(f"         first result keys: {list(res[0])[:10]}")
                sp = res[0].get("species") or {}
                if sp:
                    print(f"         species fields: {list(sp)[:8]}")
        elif isinstance(data, list) and data:
            print(f"         first item keys: {list(data[0])[:10]}")
        time.sleep(a.pause)

    print("\n3. species coverage (affects presence/absence scoring)")
    code, data, url = get("/genome/", {"per_page": 1})
    if isinstance(data, dict):
        print(f"  genomes in this release: {data.get('count')}")
    else:
        show("genome list", code, data, url)

    print("\nNext: confirm the ortholog response carries a species/taxon id "
          "per hit. That field is what the presence/absence fingerprint is "
          "built from.")


if __name__ == "__main__":
    main()
