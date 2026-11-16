#!/usr/bin/env python3
"""
Fetch protein sequences for every receptor and precursor accession.

First step toward a trained model: ESM embeddings need sequences, and this
pulls them from UniProt in batches. CPU only, no GPU, resumable.

Note on the peptide side: embed the PRECURSOR protein, not the mature peptide.
Mature peptides are 10-40 residues and protein language models are unreliable
on fragments that short -- there is too little context to condition on.

Usage:
    python fetch_sequences.py
    python fetch_sequences.py --batch 50

Outputs:
    data/sequences.fasta       all sequences
    data/sequences.csv         accession, role, length, sequence
    data/sequences_report.txt

Data: UniProtKB (CC-BY 4.0).
"""

import argparse
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

import pandas as pd

STREAM = "https://rest.uniprot.org/uniprotkb/stream"
UA = "academic script (peptide-GPCR pairing project)"


def fetch_batch(accs, timeout=180, retries=3):
    """Fetch a batch of accessions as FASTA."""
    query = " OR ".join(f"accession:{a}" for a in accs)
    params = {"query": query, "format": "fasta"}
    url = STREAM + "?" + urllib.parse.urlencode(params)
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8", "replace"), None
        except urllib.error.HTTPError as e:
            if e.code == 400:
                return None, f"HTTP 400 (bad query)"
            last = f"HTTP {e.code}"
        except Exception as e:
            last = f"{type(e).__name__}"
        time.sleep(2 * (attempt + 1))
    return None, last


def parse_fasta(text):
    """Return {accession: sequence} from UniProt FASTA headers."""
    out, acc, buf = {}, None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if acc:
                out[acc] = "".join(buf)
            # >sp|P30411|BKRB2_HUMAN ...
            parts = line[1:].split("|")
            acc = parts[1] if len(parts) > 2 else parts[0].split()[0]
            buf = []
        elif line.strip():
            buf.append(line.strip())
    if acc:
        out[acc] = "".join(buf)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default="data/pairs_resolved.csv")
    ap.add_argument("--out", default="data")
    ap.add_argument("--batch", type=int, default=50)
    ap.add_argument("--pause", type=float, default=0.5)
    a = ap.parse_args()

    log = []

    def say(m=""):
        print(m)
        log.append(m)

    df = pd.read_csv(a.pairs, dtype=str)
    role = {}
    for acc in df["receptor_uniprot"].dropna().unique():
        role[acc] = "receptor"
    for acc in df["precursor_uniprot"].dropna().unique():
        role.setdefault(acc, "precursor")
    # multi-subunit entries are pipe-joined; split and keep the parts
    expanded = {}
    for acc, r in role.items():
        for part in str(acc).split("|"):
            part = part.strip()
            if part:
                expanded.setdefault(part, r)
    role = expanded

    accs = sorted(role)
    say(f"{len(accs)} distinct accessions "
        f"({sum(1 for v in role.values() if v == 'receptor')} receptor, "
        f"{sum(1 for v in role.values() if v == 'precursor')} precursor)")

    # resume from an existing file
    path_csv = os.path.join(a.out, "sequences.csv")
    have = {}
    if os.path.exists(path_csv):
        prev = pd.read_csv(path_csv, dtype=str)
        have = dict(zip(prev["accession"], prev["sequence"]))
        say(f"resuming: {len(have)} sequences already cached")

    todo = [x for x in accs if x not in have]
    say(f"to fetch: {len(todo)}")

    failed = []
    for i in range(0, len(todo), a.batch):
        chunk = todo[i:i + a.batch]
        text, err = fetch_batch(chunk)
        if text is None:
            failed.extend(chunk)
            say(f"  batch {i // a.batch + 1}: FAILED ({err})")
            continue
        got = parse_fasta(text)
        have.update(got)
        missing = [x for x in chunk if x not in got]
        failed.extend(missing)
        say(f"  batch {i // a.batch + 1}: {len(got)}/{len(chunk)} returned"
            + (f", missing {missing}" if missing else ""))
        time.sleep(a.pause)

    if not have:
        sys.exit("no sequences retrieved")

    rows = [{"accession": acc, "role": role.get(acc, ""),
             "length": len(seq), "sequence": seq}
            for acc, seq in sorted(have.items()) if acc in role]
    out = pd.DataFrame(rows)
    os.makedirs(a.out, exist_ok=True)
    out.to_csv(path_csv, index=False)

    with open(os.path.join(a.out, "sequences.fasta"), "w") as f:
        for _, r in out.iterrows():
            f.write(f">{r.accession}|{r.role}\n")
            for j in range(0, len(r.sequence), 60):
                f.write(r.sequence[j:j + 60] + "\n")

    say()
    say(f"sequences: {len(out)}")
    for r in ("receptor", "precursor"):
        sub = out[out.role == r]
        if len(sub):
            say(f"  {r:<10} n={len(sub):>3}  length median "
                f"{int(sub.length.median())}, min {int(sub.length.min())}, "
                f"max {int(sub.length.max())}")
    short = out[out.length < 50]
    if len(short):
        say(f"  {len(short)} sequences under 50 residues — language models are "
            "unreliable on fragments this short")
    if failed:
        say(f"  unretrieved: {len(failed)} -> {failed[:10]}")

    say()
    say("Next: ESM-2 embeddings. 650M on CPU for this many proteins is slow "
        "but feasible; a GPU on Sherlock is the better route. Embed the "
        "precursor, not the mature peptide.")

    with open(os.path.join(a.out, "sequences_report.txt"), "w") as f:
        f.write("\n".join(log) + "\n")
    print(f"\nwrote {path_csv}, {a.out}/sequences.fasta and the report")


if __name__ == "__main__":
    main()
