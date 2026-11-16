# Findings

Analysis of 15 November 2026. All numbers come from `data/benchmark_report.txt`
and `data/benchmark_results.csv`, reproducible from the pipeline in the README.

## Result

Presence/absence co-evolution does not rank cognate peptide–GPCR pairs above
family-matched non-pairs.

Cold-receptor split, hard negatives, n = 1,157, base rate 0.25, three seeds.
AUROC with 95% bootstrap intervals (seed 0):

| Score | AUROC | 95% CI | AUPRC |
|---|---|---|---|
| Bit identity (Foster's metric) | 0.476 | 0.434–0.519 | 0.255 |
| Pearson correlation | 0.523 | 0.482–0.565 | 0.310 |
| Jaccard | 0.529 | 0.489–0.572 | 0.312 |
| Pearson, shared species only | 0.559 | 0.519–0.598 | 0.308 |
| **Peptide conservation alone** (`n_peptide`) | **0.578** | 0.538–0.618 | 0.324 |
| **Both-present count** (`n_both`) | **0.589** | 0.550–0.633 | 0.391 |
| **Receptor conservation alone** (`n_receptor`) | **0.629** | 0.588–0.671 | 0.466 |

Every similarity metric's interval spans 0.5. Receptor conservation — one
number per gene, carrying no information about partnership — outranks all of
them.

### Controls

- **Self-pairs: 1.0000.** Every gene scored against itself gives exactly 1.0 on
  bit identity and Jaccard, so the metric detects identity.
- **Shuffled labels: 0.501** (seeds 0.478–0.520). The evaluation does not leak.

A null result is only a measurement if the instrument works. Both controls pass.

### The signal lives in the negative set

| | Similarity metrics | Conservation confounders |
|---|---|---|
| Hard (family-matched) negatives | 0.476–0.559 | 0.578–0.629 |
| Random negatives | ~0.561 | 0.510–0.538 |

Swapping family-matched negatives for random ones makes the similarity metrics
look better and the confounders look worse. The apparent signal is a property
of which negatives are chosen, not of the pairings.

Mean bit identity is ~0.91 for cognate pairs **and** ~0.91 for hard negatives.
Foster et al. reported 0.91 for cognate pairs against 0.56 for random
protein–protein pairs. The negatives in this study reproduce their
positive-pair figure, because within the peptide–GPCR system every gene is
broadly conserved and bit identity counts shared absence as agreement.

### Why the three split schemes give identical numbers

These scores are training-free. Nothing is fitted, so no information passes
from a training fold to a test fold, and fold assignment cannot change the
answer. Cold-start splits begin to matter only once a model is trained. The
harness is infrastructure for later work, not a constraint on this result.

## How this sits in the literature

This null is the expected outcome given what is already known about
rate-covariation and profiling methods, which makes it hard to attribute to a
flawed implementation.

- Clark et al. (2012, *Genome Res* 22:714) found that direct physical
  interaction is not required to produce evolutionary rate covariation, with
  strong correlations between non-interacting but co-functional enzymes, and
  that ERC is distributed uniformly along the protein sequence rather than
  concentrated at interfaces.
- Little J, Chikina M, Clark NL (2024, *eLife* 12:RP93333) tested this directly
  on protein complexes and found no generalizable relationship between physical
  interaction and high ERC, concluding that the dominant signal comes from
  selective pressure on the protein as a whole rather than from compensatory
  change at a binding interface.

Both conclude that co-evolution measures shared **function**, not physical
**binding**. Peptide–receptor pairing is a physical binding question asked
within a set of genes that already share function — they are all parts of the
peptidergic signalling system. On that reading, the absence of discrimination
here is what the method should be expected to produce.

The contribution is that nobody had tested it for this system, and the claim in
circulation points the other way.

## What this does not show

- It does not show that peptides and their receptors lack shared evolutionary
  history. It shows that this way of measuring it carries no pairing
  information once the comparison is made within the system.
- It covers one orthology database, one pair set, and one family of metrics.
- Evolutionary rate covariation proper is untested here. The literature above
  suggests it is unlikely to separate physical partners, but that is a
  prediction, not a result. ERCnet (Forsythe et al., *Mol Biol Evol* 2025) is
  built for phylogenies with gene duplication and would be the right tool,
  given the paralogy in the chemokine families.
- A non-pair is an untested pair, not a verified negative. Only 386 of 26,082
  possible combinations are known, so some negatives are real pairings nobody
  has assayed. This caps achievable AUROC.

## A second system where the same question applies

Cho et al., *Genetic comparisons of interleukin-17 reveal a framework for
complex signaling evolution* (bioRxiv, April 2026; PMC13105070, PMID 42039560)
apply ERC to the IL-17 ligand–receptor family and report putative interactions
for IL-17D:IL-17RC and IL-17B:IL-17REL. The design has the same gap this study
was built to probe, and the evidence is in their own text.

**No statistical control is reported for the ERC analysis.** The full text
contains no occurrence of "permutation", "random", "null", "p-value" or "FDR".
A "threshold of significance" is applied to ERC scores but not derived. By
contrast, their selection analyses report explicit criteria (PAML BEB > 0.99,
MEME p < 0.05, FUBAR posterior > 0.95), so the asymmetry is specific to ERC
rather than a general looseness.

**Known pairs are evaluated by eye on an annotated heatmap.** Known
ligand–receptor pairs are marked with black boxes; unmarked high-scoring cells
are read as previously unrecognized associations. Within a family of six
ligands and six receptors every cell is already a family-matched comparison —
the regime in which the present study finds high scores regardless of
partnership.

**A canonical pair fails the test and is explained post hoc.** They report that
IL-17A showed significantly weaker ERC with its cognate receptors IL-17RA and
IL-17RC than their threshold, suggesting this reflects "functional flexibility"
arising from promiscuity. They then report IL-17E:IL-17RB as defying that
trend: also promiscuous, but with a strong ERC signal they take as evidence of
a real interaction. Weak ERC for a promiscuous ligand is explained by
promiscuity; strong ERC for a promiscuous ligand is taken as confirmation.
Both outcomes support the framework, which leaves nothing that could
contradict it.

They cite Little et al. (2024) — the paper concluding ERC predicts co-functional
rather than physical interactions — and nonetheless interpret ERC cells as
binding predictions.

**Proposed test.** Apply the evaluation in this repository to the IL-17 family:
known pairs as positives, the remaining ligand–receptor cells as
family-matched negatives, AUROC with bootstrap intervals, self-pair and
shuffled-label controls. Expected caveat: with roughly 36 cells and a handful
of positives, the intervals will be very wide, and the honest conclusion may be
that the family is underpowered for evaluation — which is itself a criticism of
drawing predictions from it.

One of their predictions has independent support: Li et al. (2026, *Nat
Immunol* 27:225–235) report IL-17REL as a decoy receptor binding IL-17A,
IL-17F and IL-17C. A decoy binding many ligands is promiscuity rather than
specific pairing, so this is weak evidence for a cognate IL-17B:IL-17REL pair
in particular.

## Secondary finding: the orthology database decides what is analysable

OMA assigns no orthology group at all (`oma_group` 0, no HOG) to 18 of 248
proteins in this dataset. All 18 are covered by Ensembl Compara, median 166
ortholog species. The gap is concentrated in the chemokine system — CCR5,
CCR2, CCL5, CCL23, CX3CL1 — rather than spread at random, and also removes
CALCR, CALCA, CALCB, CRHR1, ADCYAP1R1, AVPR2, MC1R, LGR6, HCRTR1, NPFFR2,
UTS2B, NPB and PSAP.

Any co-evolution analysis built on OMA silently excludes these genes. Evidence
is in `data/ensembl_coverage_probe.csv`.

Fingerprints from the two databases are not interchangeable: Ensembl Compara
spans a few hundred vertebrates, OMA about 3,000 species across all domains. A
cross-database comparison has to be run per database on a shared species set,
not merged.

## Data provenance note

Four defects in the GtoPdb bulk files and one of my own were found during
assembly; see the README. One is a curation error worth reporting upstream:
ligand 1866 (prokineticin-1) carries prokineticin-2's accession `Q9HC23`, while
its name, species and 86-residue sequence are all correct for prokineticin-1.
