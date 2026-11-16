# Peptide–GPCR pairing: known pairs and an evaluation harness

Public-data study on pairing peptide ligands to G protein-coupled receptors.

**Status.** Data assembly, the evaluation harness and the co-evolution
benchmark are complete. No model is fitted and nothing here predicts anything.

**Result in one line:** presence/absence co-evolution does not rank cognate
peptide–GPCR pairs above family-matched non-pairs (AUROC 0.476–0.559, all
intervals spanning chance), while receptor conservation alone scores 0.629.
Both controls pass. Full numbers, controls and interpretation in
[`FINDINGS.md`](FINDINGS.md).

## Relationship to Foster et al. (2019)

This project is downstream of Foster SR et al., *Cell* 179:895, and the overlap
is substantial enough to state up front. That paper assembled known human
peptide–GPCR pairs from GtoPdb (120 receptors, 130 ligand precursors), built
presence/absence orthology fingerprints across 313 OMA species, and reported
that cognate pairs co-evolve: average fingerprint identity of 91% for known
pairs against 56% for random protein–protein pairs.

The pair set here (386 pairs, 126 receptors, 207 peptides, 130 precursor genes)
is essentially the same set, independently assembled from a 2026 GtoPdb
release. Arriving at 130 precursor genes separately is reassurance that the
assembly is right, not a new finding.

What this project adds is an evaluation the original did not run:

- **Hard negatives.** Foster's comparison was against *random* protein–protein
  pairs, which include lineage-restricted and non-metazoan genes. Two
  vertebrate-conserved genes will look similar to each other compared with
  that, whether or not they are partners. The benchmark here uses negatives
  drawn from the right receptor family.
- **Pair-level ranking.** The fingerprint result was descriptive. Their
  predictive model was a logistic regression over 17 receptor residue
  positions, which classifies whether a receptor binds *a* peptide, not
  *which* one. Whether co-evolution discriminates between candidate partners
  is untested.
- **Scope sensitivity.** OMA now covers 2,983 species against the 313 used in
  2019. Most additions are taxa where both a GPCR and a peptide precursor are
  absent, and under bit-identity similarity shared absence counts as
  agreement. Scoring the same pairs at both scopes tests whether the reported
  signal is a conservation artifact.

The outcome may be that the signal holds, or that it collapses under hard
negatives. Both are reportable; the second is the more interesting.

## Pipeline

```
sh fetch_data.sh              # GtoPdb bulk files -> data/gtopdb/
python build_pairs_v3.py      # -> data/pairs.csv
python resolve_precursors.py  # -> data/precursors.csv   (UniProt, cross-checked)
python resolve_pairs.py       # -> data/pairs_resolved.csv
python make_splits.py --pairs data/pairs_resolved.csv --neg-ratio 3
python fetch_orthologs_all.py # OMA orthologs -> data/fingerprints_all.csv
python score_coevolution.py   # -> data/scores_*.csv
python benchmark.py           # -> data/benchmark_results.csv
```

`fetch_orthologs.py` is the same fetch restricted to 1:1 orthologs (Foster's
setting); `fetch_orthologs_all.py` drops that filter. Both are kept so the
relation type can be compared. `probe_ensembl.py` checks a second orthology
source for the proteins OMA leaves out.

`data/pairs_resolved.csv` is the canonical table. `probe_oma.py` checks the OMA
API and is a diagnostic, not part of the pipeline.

Every step prints the counts it acted on and writes a report. This is
deliberate: a bug in `resolve_pairs.py` silently deleted 99 pairs (26% of the
data) and was caught only because the pair count was printed and read. Silent
pipelines hide exactly this class of error.

## Data

IUPHAR/BPS Guide to PHARMACOLOGY, downloaded 15 November 2026. The download
page states release 2026.2; the site's current release is 2026.3.1. Record the
release when re-running. `download.jsp` requires a login, but the static files
under `/DATA/` do not; the GPCR interaction file is lowercase
(`gpcr_interactions.csv`).

Precursor accessions come from UniProtKB, resolved by exact primary gene name
on reviewed human entries.

## How a pair is defined

- **Peptide ligand**: the interaction's `Ligand ID` appears in GtoPdb
  `peptides.csv`. The `Ligand Type` column cannot be used (see defects).
- **Human receptor**: `Target Species` is exactly `Human`.
- **Human peptide**: `Ligand Species` contains `Human`; blank-species rows are
  kept and flagged (1 row).
- **Endogenous only**: `Endogenous` is `true`, so synthetic analogues are
  excluded. 477 agonist and 17 antagonist annotations survive.
  `--no-endogenous-only` relaxes this.
- **One row per (peptide, receptor)**, *not* per (precursor, receptor). One
  precursor yields several distinct peptides with different pharmacology — GCG
  gives both glucagon and GLP-1, which act at different receptors; GAST,
  CXCL12, PDYN and TAC1 likewise. The precursor accession is a feature key,
  carried alongside.

## Defects found in the source data

Each was found by a cross-check or a printed count, not by reading the files.

| Defect | Location | Effect if unhandled | Resolution |
|---|---|---|---|
| `Ligand Type` empty for every row | `gpcr_interactions.csv` | Peptide filter returns zero rows | Membership by join to `peptides.csv` |
| `Species` comma-delimited, `UniProt id` pipe-delimited | `peptides.csv` | Lists never align; fallback asserts wrong accessions | No accession asserted; resolved from UniProt |
| Gene symbols not in species order | `gpcr_interactions.csv`, ligand 1579 | `Nts` (rat) taken as the human symbol | Collapsed via accession; 1 gene affected (`P30990`) |
| Wrong accession | `peptides.csv`, ligands 1866/1867/1868 | Prokineticin-2's `Q9HC23` assigned to prokineticin-1 | Override in `overrides.csv`; reported upstream |
| Dedup keyed on precursor | `resolve_pairs.py` (ours) | 99 real pairs merged away | Keyed on peptide instead |

For the prokineticin case the name, species and 86-residue sequence are all
correct for prokineticin-1; only the accession field is wrong, which is why it
is overridden rather than excluded.

Precursor resolution agrees with a GtoPdb-listed accession for 126 of 130
genes, disagrees for 1 (the above), and has no raw value to check for 3.

### Carried-forward caveats

Some entries are not classical prohormones — HBB, PSAP, RPS19, YARS1, ANXA1,
PRNP, MT-RNR2, and the WNT and RSPO families. These are full proteins whose
fragments act at GPCRs, or Frizzled and LGR ligands. They are legitimate GtoPdb
entries but behave differently under co-evolution: a housekeeping gene is
conserved everywhere for reasons unrelated to its receptor. Results should be
stratified by classical prohormone against moonlighting protein, because a
signal arising only from the latter supports the artifact hypothesis.

## Evaluation harness

`make_splits.py` fixes what counts as a positive, a negative, and a training
row, before any feature exists.

**Negatives.** A non-pair is not a verified negative — it is an untested pair.
386 of 26,082 possible pairs are known, so 98.5% of the grid is unassayed and
some "negatives" are real pairings nobody has measured. This caps achievable
AUROC and is a property of the problem, not a bug.

- `hard`: receptor shares a family with a known partner of that peptide but is
  not itself a partner. Pool of 1,041, used in full so the negative set is
  deterministic and seed variance reflects only fold assignment.
- `random`: any non-pair. Weak; retained to quantify how much easier it is.

**Splits.**

- `cold_receptor` (primary): whole receptors held out with all their edges.
  This is the orphan condition. Folds of 215–403 rows, 67–92 positives, base
  rates 0.23–0.35.
- `cold_both` (stress test): receptors and peptides partitioned independently;
  a row is in test only if both are held out, in training only if neither is,
  and otherwise discarded. Folds of 35–76 rows — underpowered at this data
  scale, reported with that caveat rather than as a headline.
- `random_edge`: inflated comparator only. A peptide's other receptors leak
  into training.

Every scheme asserts that no held-out entity appears in training; the
assertion caught a real bug in `cold_both` during development.

Report AUPRC alongside AUROC — the base rate is about 0.25, not 0.5.

## Next

Not started, in rough order of cost:

- Apply this evaluation to the IL-17 ligand–receptor family, where an ERC
  study makes pairing predictions from an annotated heatmap without a negative
  set (see `FINDINGS.md`). Reuses the existing harness.
- Evolutionary rate covariation as a second feature, via ERCnet, which handles
  the gene duplication present in the chemokine families. The literature
  predicts it measures co-function rather than binding, so a similar null is
  likely.
- Ensembl Compara fingerprints as a parallel per-database analysis. Not
  mergeable with the OMA matrix — different species scope.
- Receptor and precursor sequence embeddings with a cold-start
  link-prediction model, for which the split harness already exists.

## Licence and credit

GtoPdb data are licensed ODbL, contents CC-BY-SA 4.0.
Harding SD, Armstrong JF, Faccenda E, Southan C, Alexander SPH, Davenport AP,
Spedding M, Davies JA (2026) The IUPHAR/BPS Guide to PHARMACOLOGY in 2026.
*Nucleic Acids Res* 54(D1):D1446–D1456. doi:10.1093/nar/gkaf1067

UniProtKB: The UniProt Consortium. CC-BY 4.0.

OMA: Altenhoff A et al. (2024) *Nucleic Acids Res* 52(D1):D513–D521.

Study design after Foster SR et al. (2019) *Cell* 179:895.
Baseline co-evolution tool: Ji F et al. (2022) *Nucleic Acids Res*
50(W1):W246–W253 (DEPCOD).
