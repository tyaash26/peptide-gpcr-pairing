#!/usr/bin/env python3
"""
Build a self-contained, searchable HTML page of the peptide-GPCR pairs.

Reads pairs_resolved.csv and (if present) the scored pairs, embeds them as JSON
in a single HTML file, and writes it out. No external requests at runtime, so
the page works offline and can be published or hosted anywhere.

Usage:
    python build_search_page.py
    python build_search_page.py --scores data/class3_scored --out pairs_browser.html

Data: IUPHAR/BPS Guide to PHARMACOLOGY (ODbL / CC-BY-SA 4.0) and OMA.
"""

import argparse
import glob
import json
import os

import pandas as pd

TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<title>Peptide–GPCR pairs: browser</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Spectral:ital,wght@0,300;0,400;0,600&display=swap" rel="stylesheet">
<style>
  :root{
    box-sizing:border-box;
    padding-top:env(safe-area-inset-top,0px);
    padding-bottom:env(safe-area-inset-bottom,0px);
    --paper:#fff; --ink:#16222b; --soft:#56646e; --rule:#dde2e5;
    --wash:#f4f6f6; --teal:#2e6e6b; --flag:#b0442e; --hit:#e8f2f1;
  }
  @media (prefers-color-scheme: dark){
    :root:not([data-theme="light"]){
      --paper:#0f161b; --ink:#e6eaec; --soft:#9aa7af; --rule:#28333b;
      --wash:#172128; --teal:#6fb3ad; --flag:#e08466; --hit:#1d2e2c;
    }
  }
  :root[data-theme="dark"]{
    --paper:#0f161b; --ink:#e6eaec; --soft:#9aa7af; --rule:#28333b;
    --wash:#172128; --teal:#6fb3ad; --flag:#e08466; --hit:#1d2e2c;
  }
  *,*::before,*::after{box-sizing:inherit}
  html,body{height:100%}
  body{margin:0;background:var(--paper);color:var(--ink);
    font-family:Spectral,Georgia,serif;font-weight:300;font-size:16px;
    line-height:1.5;font-variant-numeric:tabular-nums lining-nums}
  .wrap{max-width:62rem;margin:0 auto;padding:0 1.2rem 5rem}
  header{padding:2.5rem 0 0}
  h1{font-size:1.6rem;font-weight:600;margin:0 0 .3rem;letter-spacing:-.01em}
  .sub{color:var(--soft);font-style:italic;margin:0 0 1.5rem}
  .controls{position:sticky;top:env(safe-area-inset-top,0px);
    background:var(--paper);padding:.8rem 0;border-bottom:1px solid var(--rule);
    z-index:5}
  input[type=search],select{font-family:inherit;font-size:1rem;
    padding:.5rem .7rem;border:1px solid var(--rule);background:var(--wash);
    color:var(--ink);border-radius:2px}
  input[type=search]{width:100%;max-width:26rem}
  select{margin-top:.5rem}
  .row2{display:flex;gap:.6rem;flex-wrap:wrap;align-items:center;
    margin-top:.5rem}
  .count{color:var(--soft);font-size:.9rem}
  .scroll{overflow-x:auto;margin-top:.4rem}
  table{border-collapse:collapse;width:100%;font-size:.9rem}
  th,td{text-align:left;padding:.4rem .7rem .4rem 0;
    border-bottom:1px solid var(--rule);white-space:nowrap}
  th{font-weight:600;cursor:pointer;user-select:none;position:relative}
  th:hover{color:var(--teal)}
  th[aria-sort]{color:var(--teal)}
  td.wrapcell{white-space:normal;min-width:14rem}
  mark{background:var(--hit);color:inherit;padding:0 .1em}
  .empty{padding:2.5rem 0;color:var(--soft);font-style:italic}
  footer{margin-top:3rem;padding-top:1rem;border-top:1px solid var(--rule);
    font-size:.82rem;color:var(--soft)}
  .pos{color:var(--teal)}
  .neg{color:var(--soft)}
  @media(max-width:40rem){th,td{padding-right:.5rem;font-size:.85rem}}
</style>
</head>
<body>
<div class="wrap">
<header>
  <h1>Peptide–GPCR pairs</h1>
  <p class="sub">__NPAIRS__ curated human pairs · __NREC__ receptors ·
  __NPEP__ peptides. Search any field; click a column to sort.</p>
</header>

<div class="controls">
  <input type="search" id="q" placeholder="Search receptor, peptide, gene, family, accession…" autocomplete="off">
  <div class="row2">
    <select id="fam"><option value="">All receptor families</option></select>
    <select id="act"><option value="">All actions</option></select>
    <span class="count" id="count"></span>
  </div>
</div>

<div class="scroll"><table id="tbl">
  <thead><tr id="hrow"></tr></thead>
  <tbody id="body"></tbody>
</table></div>
<div class="empty" id="empty" hidden>Nothing matches that search.</div>

<footer>
  <p>Pair data: IUPHAR/BPS Guide to PHARMACOLOGY, ODbL; contents CC-BY-SA 4.0.
  Harding SD et al. (2026) <i>Nucleic Acids Res</i> 54(D1):D1446–D1456.
  Orthology counts from OMA, Altenhoff A et al. (2024).</p>
  <p id="scorenote"></p>
</footer>
</div>

<script>
const DATA = __DATA__;
const COLS = __COLS__;
const HAS_SCORES = __HASSCORES__;

const body = document.getElementById('body');
const hrow = document.getElementById('hrow');
const q = document.getElementById('q');
const famSel = document.getElementById('fam');
const actSel = document.getElementById('act');
const countEl = document.getElementById('count');
const emptyEl = document.getElementById('empty');

let sortKey = null, sortDir = 1;

COLS.forEach(c => {
  const th = document.createElement('th');
  th.textContent = c.label;
  th.dataset.key = c.key;
  th.addEventListener('click', () => {
    if (sortKey === c.key) { sortDir = -sortDir; }
    else { sortKey = c.key; sortDir = 1; }
    [...hrow.children].forEach(x => x.removeAttribute('aria-sort'));
    th.setAttribute('aria-sort', sortDir === 1 ? 'ascending' : 'descending');
    render();
  });
  hrow.appendChild(th);
});

function uniqueSorted(key) {
  return [...new Set(DATA.map(r => r[key]).filter(v => v != null && v !== ''))]
    .sort((a, b) => String(a).localeCompare(String(b)));
}
uniqueSorted('receptor_family').forEach(v => {
  const o = document.createElement('option'); o.value = v; o.textContent = v;
  famSel.appendChild(o);
});
uniqueSorted('action').forEach(v => {
  const o = document.createElement('option'); o.value = v; o.textContent = v;
  actSel.appendChild(o);
});

function esc(s) {
  return String(s).replace(/[&<>"]/g, c =>
    ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
}
function highlight(text, needle) {
  const t = esc(text == null ? '' : text);
  if (!needle) return t;
  const i = t.toLowerCase().indexOf(needle.toLowerCase());
  if (i < 0) return t;
  return t.slice(0, i) + '<mark>' + t.slice(i, i + needle.length) +
         '</mark>' + t.slice(i + needle.length);
}

function matches(r, needle) {
  if (!needle) return true;
  const n = needle.toLowerCase();
  return COLS.some(c => String(r[c.key] ?? '').toLowerCase().includes(n));
}

function render() {
  const needle = q.value.trim();
  const fam = famSel.value, act = actSel.value;
  let rows = DATA.filter(r =>
    matches(r, needle) &&
    (!fam || r.receptor_family === fam) &&
    (!act || r.action === act));

  if (sortKey) {
    const col = COLS.find(c => c.key === sortKey);
    rows = rows.slice().sort((a, b) => {
      let x = a[sortKey], y = b[sortKey];
      if (col && col.num) {
        x = (x === '' || x == null) ? -Infinity : Number(x);
        y = (y === '' || y == null) ? -Infinity : Number(y);
        return (x - y) * sortDir;
      }
      return String(x ?? '').localeCompare(String(y ?? '')) * sortDir;
    });
  }

  countEl.textContent = rows.length === DATA.length
    ? DATA.length + ' pairs'
    : rows.length + ' of ' + DATA.length + ' pairs';

  emptyEl.hidden = rows.length > 0;
  body.innerHTML = rows.map(r => '<tr>' + COLS.map(c => {
    const v = r[c.key];
    const cls = c.wrap ? ' class="wrapcell"' : '';
    return '<td' + cls + '>' + highlight(v, needle) + '</td>';
  }).join('') + '</tr>').join('');
}

q.addEventListener('input', render);
famSel.addEventListener('change', render);
actSel.addEventListener('change', render);
render();

if (HAS_SCORES) {
  document.getElementById('scorenote').textContent =
    'Co-evolution columns are presence/absence fingerprint scores over ' +
    'OMA orthologs. They do not separate true pairs from matched non-pairs ' +
    'strongly: AUROC about 0.56 under the primary design. Treat a high score ' +
    'as weak evidence at best.';
}
</script>
</body>
</html>
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default="data/pairs_resolved.csv")
    ap.add_argument("--scores", default="data/class3_scored")
    ap.add_argument("--out", default="pairs_browser.html")
    a = ap.parse_args()

    df = pd.read_csv(a.pairs, dtype=str).fillna("")
    print(f"{len(df)} pairs from {a.pairs}")

    cols = [
        {"key": "receptor_gene", "label": "Receptor"},
        {"key": "receptor_family", "label": "Family", "wrap": True},
        {"key": "peptide_name", "label": "Peptide", "wrap": True},
        {"key": "precursor_gene", "label": "Precursor"},
        {"key": "action", "label": "Action"},
        {"key": "receptor_uniprot", "label": "Receptor UniProt"},
        {"key": "precursor_uniprot", "label": "Precursor UniProt"},
    ]

    # attach co-evolution scores for the true pairs, if available
    has_scores = False
    score_file = sorted(glob.glob(os.path.join(a.scores, "scores_*seed0.csv")))
    if score_file:
        s = pd.read_csv(score_file[0])
        s = s[s["label"] == 1]
        keep = {"receptor_gtopdb_id": "receptor_gtopdb_id",
                "peptide_gtopdb_id": "peptide_gtopdb_id",
                "all__bit_identity": "bit_identity",
                "shared__jaccard": "jaccard",
                "all__n_receptor": "n_rec_orthologs",
                "all__n_peptide": "n_pep_orthologs"}
        s = s[[c for c in keep if c in s.columns]].rename(columns=keep)
        for c in ("bit_identity", "jaccard"):
            if c in s.columns:
                s[c] = s[c].astype(float).round(3)
        df["receptor_gtopdb_id"] = df["receptor_gtopdb_id"].astype(str)
        df["peptide_gtopdb_id"] = df["peptide_gtopdb_id"].astype(str)
        s["receptor_gtopdb_id"] = s["receptor_gtopdb_id"].astype(str)
        s["peptide_gtopdb_id"] = s["peptide_gtopdb_id"].astype(str)
        df = df.merge(s, on=["receptor_gtopdb_id", "peptide_gtopdb_id"],
                      how="left")
        n_scored = df["bit_identity"].notna().sum() if "bit_identity" in df else 0
        print(f"  attached co-evolution scores to {n_scored} pairs "
              f"(from {os.path.basename(score_file[0])})")
        if n_scored:
            has_scores = True
            cols += [
                {"key": "bit_identity", "label": "Bit identity", "num": True},
                {"key": "jaccard", "label": "Jaccard", "num": True},
                {"key": "n_rec_orthologs", "label": "Rec. orthologs",
                 "num": True},
                {"key": "n_pep_orthologs", "label": "Pep. orthologs",
                 "num": True},
            ]
    else:
        print(f"  no score files under {a.scores}; building without scores")

    keys = [c["key"] for c in cols]
    records = df[[k for k in keys if k in df.columns]].astype(object).where(df[[k for k in keys if k in df.columns]].notna(), "")
    records = records.astype(str).replace({"nan": ""}).to_dict("records")

    html = (TEMPLATE
            .replace("__DATA__", json.dumps(records, separators=(",", ":")))
            .replace("__COLS__", json.dumps(cols))
            .replace("__HASSCORES__", "true" if has_scores else "false")
            .replace("__NPAIRS__", f"{len(df):,}")
            .replace("__NREC__", str(df["receptor_gene"].nunique()))
            .replace("__NPEP__", str(df["peptide_name"].nunique())))

    with open(a.out, "w", encoding="utf-8") as f:
        f.write(html)
    size = os.path.getsize(a.out) / 1024
    print(f"wrote {a.out} ({size:.0f} KB, self-contained)")


if __name__ == "__main__":
    main()
