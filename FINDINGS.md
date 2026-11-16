# Findings

Analysis of 15 November 2026. All numbers come from `data/benchmark_report.txt`
and `data/benchmark_results.csv`, reproducible from the pipeline in the README.

## Result

**The answer depends on how negatives are constructed, and that dependence is
the finding.** Under same-family negatives, presence/absence co-evolution does
not discriminate cognate pairs at all. Under degree-matched negatives drawn
from the wider peptide-receptor class, a small but statistically robust effect
appears (AUROC ≈ 0.56, permutation p < 0.002 across three seeds). Both results
come from the same pairs, the same fingerprints and the same metrics.

| Negative set | n | Bit identity AUROC | Verdict |
|---|---|---|---|
| Same family, unmatched degree | 1,157 | 0.476 [0.434–0.519] | no discrimination; degree leaks (0.690) |
| Same family, degree-matched | 520 | 0.434 | unusable — 252 of 386 positives have no match |
| Class-wide, degree-matched 1:1 | 567 | 0.563, p = 0.018 | borderline; one metric of six |
| **Class-wide, degree-matched 3:1** | **1,104** | **0.573, p = 0.0005** | **clears correction; metrics agree** |

The 3:1 class-matched design is the only one that is both adequately powered
and free of the degree leak, so it is the primary analysis. The sections below
were written against the earlier same-family result and are retained because
the reasoning still applies to that design.

### Primary analysis (class-wide, degree-matched, 3:1)

Pre-registered before the numbers were seen: primary metric
`all__bit_identity`, primary statistic the permutation p-value, significance
threshold p < 0.0083 (Bonferroni for six metrics).

Permutation null built by shuffling labels within folds, 2,000 permutations,
three seeds (`permutation_test.py`):

| Metric | AUROC (3 seeds) | p (3 seeds) | Clears correction |
|---|---|---|---|
| Bit identity | 0.573, 0.559, 0.559 | 0.0005, 0.0015, 0.0015 | yes, all |
| Pearson | 0.563, 0.553, 0.555 | 0.0025, 0.0035, 0.0025 | yes, all |
| Jaccard | 0.560, 0.550, 0.551 | 0.0040, 0.0080, 0.0070 | yes, all |
| Pearson (shared) | 0.538, 0.541, 0.542 | 0.079, 0.018, 0.025 | no |
| Receptor conservation | 0.494, 0.501, 0.497 | 0.75, 0.94, 0.80 | no (confounder, as intended) |
| Peptide conservation | 0.498, 0.506, 0.506 | 0.97, 0.73, 0.77 | no (confounder, as intended) |

Null centres sit at 0.496–0.503, so the theoretical chance line was correct
here; an earlier single-shuffle control returning 0.458 was noise at smaller n.

Three metrics clearing correction in all three seeds, with both conservation
confounders neutralised, is the pattern a real effect produces. One metric of
six at p ≈ 0.02 — which is what the 1:1 design gave — is not.

### Not a family-separation artifact

The recurring confound in this project is a metric separating receptor families
rather than identifying pairs. Per-family AUROC for bit identity, families with
at least 40 rows:

opioid 0.778 · tachykinin 0.795 · somatostatin 0.573 · chemokine 0.551 ·
neuropeptide Y 0.532 · cholecystokinin 0.418 · relaxin 0.398 · formylpeptide
0.389

Heterogeneous in both directions, and no single family carries the pooled
value, so the effect is not produced by between-family structure. It also is
not uniform: three of eight families run the other way.

### What the effect is worth

AUROC 0.56 is statistically solid at n ≈ 1,100 and practically small. It does
not support ranking candidate ligands for an orphan receptor, and no claim of
that kind should be made from it. The defensible statement is that a weak
signal exists under this design and is absent under a stricter one.

### Superseded framing

An earlier version of this document concluded flatly that co-evolution does not
discriminate cognate pairs. That conclusion holds for same-family negatives and
is wrong as a general claim. The change came from fixing the degree leak and
adding power, not from selecting a favourable metric: the metric, statistic and
threshold were fixed in advance of the 3:1 run.

## Earlier analysis: same-family negatives

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

**On the 0.629: it is not a biological confounder.** An earlier version of this
document read it as evidence that the similarity metrics track gene
conservation. That interpretation was tested three ways and is wrong; see
"Taking the conservation effect apart" below. In the primary (3:1
degree-matched) analysis the conservation confounders sit at chance.

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

## Taking the conservation effect apart

Receptor conservation at AUROC 0.629 looked like a confounder the similarity
metrics were tracking. Three explanations were tested; all three fail.
Reproduce with `benchmark_strata.py`.

**Not driven by moonlighting proteins.** Nineteen precursor genes are not
classical prohormones (housekeeping and structural proteins whose fragments act
at GPCRs; the set is hand-curated in `analyze_network.py`). Excluding them
moves receptor conservation from 0.629 to 0.636. They are 35 of 386 pairs, 9%
of the network.

**Not a proxy for receptor promiscuity.** Spearman correlation between a
receptor's ortholog count and its number of known peptide partners is 0.064.

**It is a pooling artifact.** Computed within a single receptor family the
effect disappears: AUROC 0.423 among chemokine receptors (n = 798) and 0.480
among the glucagon family (n = 54), against 0.629 pooled. Only those two
families have enough rows to evaluate. Families whose receptors have more
orthologs carry a higher positive rate in the negative set, so pooling across
families manufactures a correlation present in none of them — Simpson's
paradox. Conservation carries no within-family pairing information either.

## A flaw in family-matched negative sampling

Scoring candidate pairs by receptor degree alone — how many peptides that
receptor is already known to bind, which says nothing about any particular
peptide — gives AUROC 0.690, above receptor conservation (0.629) and above
every similarity metric (0.476–0.559).

Negatives are drawn as same-family receptors a peptide does not bind. A
receptor with many known partners is therefore more likely to be the positive
for any peptide in its family, so degree predicts the label by construction.
Any benchmark sampling negatives per-entity from a curated interaction set
inherits this. It is independent of the conservation effect (Spearman 0.064).

This does not weaken the main result — the similarity metrics score *below*
both label-correlated baselines — but it means the absolute AUROC values here
are not a clean measure of anything, and a future design should either
condition on degree or sample negatives to balance it.

### Where the sampling breaks down entirely

For small, exclusive families it cannot generate negatives. Among WNT and RSPO
ligands, Frizzled and LGR receptors appear as 11 positives and zero negatives,
because nearly every family member is a true partner. A post-hoc subgroup
analysis of the moonlighting stratum returned AUROC 0.813 [0.683–0.917] for
Pearson correlation, which on inspection is a classifier separating
Frizzled/LGR genes from chemokine genes rather than detecting pairs: that
stratum contains 18 chemokine negatives against 1 positive, and 11 Frizzled/LGR
positives against 0 negatives. Per-family results are uninterpretable wherever
most family members are genuine partners.

## Hypotheses tested and discarded

Recorded because the main result's credibility rests on them having been
tested rather than assumed.

| Hypothesis | Test | Outcome |
|---|---|---|
| OMA's `rel_type=1:1` filter removes paralog-rich families | Refetched without the filter | Wrong — the same 18 proteins have no orthology group either way |
| Moonlighting proteins drive the conservation effect | Stratified benchmark | Wrong — 0.629 vs 0.636 excluded |
| Conservation is a proxy for receptor promiscuity | Spearman vs degree | Wrong — 0.064 |
| The moonlighting stratum shows real co-evolution signal | Family composition of those rows | Wrong — family separation, not pairing |
| The null holds under any negative-set design | Degree-matched class-wide negatives at 3:1 | Wrong — a weak effect (AUROC 0.56, p < 0.002) appears once the degree leak is removed and power is adequate |

## Network structure

From `analyze_network.py`. 386 pairs, 126 receptors, 207 peptides, 129
precursor genes.

- **Peptides are specific, receptors promiscuous.** 53% of peptides bind
  exactly one receptor (median 1, mean 1.86, max 7). Only 24% of receptors have
  a single peptide (median 2, mean 3.06, max 13, CCR3).
- **Chemokines dominate.** 21 receptors, 46 peptides, 88 pairs — and 798 of the
  1,157 scored rows. This is substantially a chemokine benchmark, and
  chemokines are the family OMA cannot fingerprint.
- **32 of 129 precursors yield several distinct peptides** (GAST 7, CXCL12 6,
  PDYN 6, GCG 5), which is why pairs are keyed on the peptide.
- Most pairs per receptor by family: opioid 6.67, calcitonin 5.0, chemokine
  4.19.

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
