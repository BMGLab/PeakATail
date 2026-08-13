# Laughney LUAD — 3′UTR switch genes: function + miRNA mechanism

From OUR real PeakATail switch calls (`switches_annotated.csv`, atlas + no-atlas),
ENSG→symbol via the local GRCh38.99 GTF (exact, not guessed). Two questions:
(1) which genes really switch/shorten/lengthen and what do they do; (2) does the
**lost distal 3′UTR** carry miRNA sites — the de-repression mechanism.

> **Read the caveats first — they bound every claim.**

---

## 0. Status after the corrected sweep (August 2026) — READ THIS FIRST

This document was written against the **pre-fix** switch calls. The corrected
cohort sweep (`RERUN_2026-08_fixed`) has since tested the two load-bearing
assumptions underneath it, and **both fail**. The gene-level narrative below is
retained as a record of the hypothesis-generation exercise, but the following
supersedes it. See `technical_report.md` §12 for the full analysis and figures.

**(1) The stage direction is not real.** Mean PDUI per (cell type, stage)
correlates with the fraction of zero-coverage (gene, cell) rows at
**Pearson r = −0.997** (p = 3.9 × 10⁻⁸³, n = 76). Conditioning on ≥1 read flattens
the Normal → Met trajectory (0.387 / 0.392 / 0.370 / 0.380); at ≥5 reads it
*reverses* (0.416 / 0.430 / 0.422 / 0.426) and only 4 of 18 cell types still trend
down, versus 15 of 18 unconditionally. **Every "shortening" and "lengthening" call
in the tables below inherits this.** Direction was the one thing the original
document said to trust; it is not trustworthy.

**(2) The miRNA/ARE de-repression mechanism is not supported.** Scanning the
actual lost distal segment for all 3,641 genes with a proximal/distal PAS pair,
against a length-matched proximal control (`technical_report.md` §12.5):

- miRNA 7mer-m8 seed density is **0.510/kb lost vs 0.520/kb control** — no
  enrichment (Wilcoxon p = 0.058), and the sign points the wrong way.
- The ARE excess (1.70 vs 1.60 pentamers/kb) is ~6% relative and fully explained
  by AT content (ρ = +0.82 with AT fraction); the functional `WTTTATTTAW` nonamer
  shows no difference (p = 0.63).
- **The interval is mostly not 3′UTR**: median proximal-to-distal span is
  **10,779 bp**, 80% exceed 3 kb. The classic PDUI layer is pairing PAS too far
  apart to be a tandem 3′UTR pair, so §2's gene-by-gene reasoning about "the lost
  distal 3′UTR" is reasoning about intervals that are largely not UTR.

**(3) Significance counts index power, not biology.** 43.9% of all 1,233,328
Fisher tests are "significant" at FDR < 0.05 — the pseudoreplication signature.
Use hit counts to rank, never to assert a per-gene result.

**What survives.** The recurrence structure is still informative as a *candidate
list*, now with symbols resolved from the run's own GTF: 8,441 genes carry a trend
call, 3,007 (35.6%) in ≤2 cell types, and the most recurrent are **EZR, NCOA3,
CEP170, USP12** (19 of ≤24 cell types), then **NMT1, MEF2A, MAP4K3, RP2, MZT1,
RILPL2, ACOT9, PTPRJ** (18). But recurrence tracks **testability**, not shared
regulation — private genes actually have *larger* mean |slope| than recurrent ones
(0.0137 vs 0.0104, Mann–Whitney p = 0.029). Note that **MEF2A** appears in the
recurrent set here and in the original "not in the Cancer Gene Census" list below.

**What it would take to revive a biological claim** (in order):
1. Restrict PDUI to proximal/distal pairs inside the same annotated 3′UTR.
2. Test per gene only where both compared stages have adequate coverage, with
   depth as an explicit covariate and library-size-matched cells.
3. Re-derive direction on that subset before any gene is named in a figure.


## Caveats (project-verified, non-negotiable)
1. **No gene passes the stagewise test.** `differs_across_stages=False` for 100% of
   rows; `omnibus_q` never significant. These are **descriptive, hypothesis-
   generating** switches, not validated stage effects.
2. **PDUI saturates 0/1** and top calls are non-monotonic single-arm flips — the
   detection-rate / "98% PDUI-zeros" artifact. Trust **direction**, distrust
   magnitude; prefer the few **sub-saturation, multi-stage-graded** calls.
3. **Atlas precision is circular** → the **no-atlas table is the honest one**;
   genes in BOTH tables are most trustworthy.
4. **Compartment**: most switches are in **immune/stromal** cells (macrophage_M2,
   cytotoxic_T/NK, DC, ciliated), NOT tumor epithelium. Only **AE2 / EPITHELIAL**
   calls are "tumor-cell APA."
5. Direction inherits **fisher pseudoreplication over reads**; subset cohort;
   signature-scored cell types.

---

## 1. Curated switch genes with function (credibility × cancer relevance)

### Shortening (proximal PAS gained; PDUI drops Normal→met) — canonical cancer direction
| Gene | Compartment | ΔPDUI | Function | Coherence |
|---|---|---|---|---|
| **MDM2** | cytotoxic_T/NK (+6) | −1.0 | E3 ligase, master **p53** negative regulator | **Coherent** — de-represses an oncogene |
| **HK1** | DC/endothelial (7 CT) | −1.0 | Hexokinase-1, glycolysis (Warburg) | Coherent, broad support |
| **GLS** | DC/macrophage | −1.0 | Glutaminase, glutaminolysis (druggable, CB-839) | Coherent (NSCLC glutamine addiction) |
| **FNDC3B** | ciliated | −0.98 | EMT/invasion oncogene | Coherent |
| **DIAPH1** | macrophage_M2 | −0.84 | Formin, motility | Coherent; **sub-saturation → cleaner** |
| **AKAP13** | ciliated | −1.0 | Rho-GEF oncogenic scaffold | Coherent |
| **NUDT4** | macrophage_M2 + **epithelial** | −1.0 | Nudix phosphohydrolase | **Atlas+noatlas** (high confidence); direction neutral |
| ACADM | DC + **epithelial** | −1.0 | Fatty-acid β-oxidation | Atlas+noatlas; metabolic |

### Lengthening (distal PAS gained; against the cancer grain)
| Gene | Compartment | ΔPDUI | Function | Coherence |
|---|---|---|---|---|
| **MSI2** | macrophage_M2 + **epithelial** | +1.0 (brain-met) | RNA-BP; stabilizes oncogenic mRNAs; NSCLC met driver | Surprising (onco lengthening) but MSI2 is itself a UTR regulator |
| **CD44** | cytotoxic_T/NK | +1.0 (bone-met) | CSC/metastasis marker, miRNA ceRNA hub | Surprising; heavily miRNA-regulated UTR |
| **LATS2** | monocytic/T | +1.0 | Hippo kinase, **tumor suppressor** | **Coherent (pro-tumor)** — repress a TS |
| **ARID1B** | exhausted_reg_T | +0.70 | SWI/SNF, **tumor suppressor** | Coherent; **sub-saturation → cleaner** |
| **CAMTA1** | DC/endothelial | +0.93 | 1p36 **tumor suppressor** | Coherent; **7-stage graded → cleanest** |
| **SH3GLB1** | monocytic | +0.83 | Bif-1, apoptosis TS | Coherent |

**Coherent set** (fits the model): shortening of oncogenes/metabolic-motility
effectors (MDM2, HK1, GLS, FNDC3B, DIAPH1, AKAP13); lengthening of tumor
suppressors (LATS2, ARID1B, CAMTA1, SH3GLB1). **Against-grain** (scrutinize):
oncogene lengthening (MSI2, CD44, REL, YTHDF3) — several are RNA-regulatory genes
(plausible autoregulatory UTR remodeling).

Cleanest (sub-saturation, multi-stage): **DIAPH1, ARID1B, CAMTA1, GLS, ANXA11, SH3GLB1.**

---

## 2. Does the lost distal 3′UTR carry miRNA sites? (the mechanism)

Mechanism (Mayr & Bartel 2009; Sandberg 2008): shortening removes distal miRNA/ARE
sites → mRNA de-repressed → more protein. Checked our shortening genes' 3′UTRs
against TargetScan / miRTarBase / luciferase-validated literature. **Honest limit:
web tools can't resolve whether a validated site sits inside OUR exact interpeak
interval — reasoning is at 3′UTR level.**

**Tier A — consistent with oncogenic de-repression:**
- **CD44 (AE2)** — the single cleanest candidate. Oncogene; 3′UTR carries
  **luciferase-validated miR-34a (2 sites), miR-373-3p, miR-328-3p**; in the
  tumor-relevant AE2 compartment. *Undercut* by CD44 also being called
  *lengthening* in ciliated cells (saturated ±1 → direction not robust).
- **MECOM/EVI1 (macrophage)** — oncogene, **validated miR-22 site** — but wrong
  compartment + leukemic context.

**Tier B — the important honest finding: the best-validated miRNA biology in our
shortening list is on TUMOR SUPPRESSORS**, so de-repression would be *anti*-oncogenic:
- **CDKN1B/p27 (AE2)** — textbook **miR-221/222** double site (validated). Shortening → ↑p27 → arrest (opposite sign).
- **NF1 (macrophage)** — RAS-GAP TS; validated **miR-193b, miR-103, miR-107**. ↑NF1 → less RAS (opposite sign).
- **TIMP2 (NK/T)** — metastasis suppressor; **miR-197-3p luciferase-validated *specifically in LUAD metastasis*** (the most disease-matched miRNA in the panel) + miR-210. ↑TIMP2 → anti-metastatic (opposite sign).
- **TGFBR2 (ciliated)** — validated miR-301a/miR-93/miR-135b; sign genuinely ambiguous (TGF-β dual role).

**Lengthening genes behave as predicted** (gain distal region → gain repression):
NR3C1 (validated miR-124/miR-142-3p), ATG7 (validated miR-17/miR-375) — but every
top ciliated lengthening call shows the identical `Normal 0 → MetBone 1` flip = an
artifact signature, not independent biology.

---

## 2.5 Isoform switches (protein-changing) — beyond 3′UTR length

Not every switch is 3′UTR-only. Where the gene switches which **transcript/isoform**
is dominant (`utr_class=multi_isoform`, `n_utr_isoforms≥2`), the *protein* can
change. High-confidence core = **atlas ∩ no-atlas = 61 genes**; classified against
Ensembl/APPRIS/UniProt/Pfam:

| Gene (compartment) | Switch class | Real protein change? | Consequence |
|---|---|---|---|
| **BMPR2** (immune-biased) | **coding — C-term truncation** | **Yes** | short isoform keeps kinase domain but **deletes the ~508-aa cytoplasmic regulatory tail** |
| **NF1** (immune/epithelial) | **coding — in-frame exon-23a cassette** | **Yes** | 21-aa cassette **inside the RAS-GAP domain** — tunes GAP catalysis (on-mechanism for the TSG) |
| **YAP1** (immune) | **coding — WW-domain #** | **Yes** | 2-WW vs 1-WW isoform → different PPxY partner binding (textbook YAP switch) |
| **CD46** (epithelial/immune) | **coding — alt last exon** | **Yes** | mutually-exclusive cytoplasmic tails CYT-1 vs CYT-2 (different signalling) |
| **CYLD** (immune) | **3′UTR-only** | No (regulatory) | byte-identical 953-aa protein — the clean "same protein, different regulation" case |
| TPM4 / CREB1 | mixed / NMD-contaminated | doubtful | 2/3 flagged transcripts are NMD biotype |

**Methodological caveat surfaced (important):** **25% of flagged "dominant
transcripts" are Ensembl `nonsense_mediated_decay` biotype**, enriched for
RNA-binding/splicing factors (HNRNPU ×7, LUC7L3 ×5, RSRP1 ×4, CREM, PDIA3, FUS,
HNRNPC…). Either genuine AS-NMD autoregulation, or PeakATail matching a 3′ peak to
an NMD isoform's terminal exon — **audit before any isoform claim**. Also: no
`length_periso` per-transcript files exist in the report dirs, so per-isoform
length can't be quantified from this run; and these are cell-type/stage dominance
differences under saturated PDUI, not proven tumour-clonal isoform switches.

## 2.6 Database validation (COSMIC CGC / OncoKB / miRTarBase / TargetScan)

**Cancer-gene status — validated:**
- **Bona-fide drivers** (COSMIC Tier-1 + OncoKB): **NF1** (TSG, mutated ~8–11% LUAD), **MDM2** (oncogene), **YAP1** (oncogene), **CYLD/CDKN1B/LATS2/ARID1B** (TSGs), **MSI2** (oncogene).
- **NOT in the Cancer Gene Census** (do **not** call these "drivers"): **CD44, BMPR2, MEF2A, HK1, GLS, TIMP2** — CD44 is a stem/metastasis *marker*, HK1/GLS metabolic, TIMP2 an ECM regulator.

**miRNA–3′UTR interactions — validated vs contested:**
| Pair | Verdict | Note |
|---|---|---|
| **CDKN1B / miR-221-222** | **Solid** (luciferase dual-mutant + WB + in vivo) | gold-standard |
| **TIMP2 / miR-197-3p** | **Solid — and validated in LUAD** | the most disease-matched interaction |
| **NF1 / miR-193b** | Validated (HNSCC context) | tissue extrapolation to LUAD |
| **CD44 / miR-34a** | **CONTESTED** ⚠ | independent **replication failed** (eLife 2018); TargetScan site poorly conserved — **do not present as settled** |
| **MECOM / miR-22** | caveated ⚠ | founding paper **retracted**; cite the 2025 re-validation instead |

> This **corrects the earlier "CD44 is the cleanest candidate" framing**: CD44's
> miR-34a link failed replication. The genuinely solid, disease-matched miRNA
> mechanism in our set is **TIMP2/miR-197-3p** (and CDKN1B/miR-221-222) — but both
> are tumour-suppressor-direction, so shortening them would be *anti*-oncogenic.

## 3. Bottom line
- The mechanism is **real and well-precedented**, but **our data cleanly supports
  it for only one gene — CD44 (AE2)** (oncogene ↑ via loss of validated distal
  miR-34a/miR-373) — and even that is undercut by a contradictory call in another
  cell type.
- **A blanket "3′UTR shortening → oncoprotein de-repression" claim is NOT
  supported** by our gene set: the shortening genes with the best-validated
  distal-miRNA biology are **tumor suppressors** (p27, NF1, TIMP2), where
  de-repression opposes tumorigenesis.
- Two structural caveats dominate: most switches are **immune/stromal, not tumor
  epithelium**, and **PDUI is saturated 0/1 with pseudoreplicated direction** —
  consistent with the project's prior finding that "global 3′UTR shortening" here
  is largely a detection-rate artifact.

## 4. What to do before putting a mechanism in a figure
1. **Re-call in the epithelial/AE2 compartment on non-saturated PDUI** (fix the
   detection-rate artifact) — only then are per-gene directions defensible.
2. **Map exact miRNA-site coordinates against the interpeak interval** (local
   TargetScan/scan of the distal segment, not web lookup) to prove the site is
   actually lost.
3. **Highest-value candidates**: **CD44** (oncogene, tumor compartment),
   **TIMP2** (miR-197-3p LUAD-validated), **CDKN1B** (miR-221/222) — the three
   with disease-matched, luciferase-validated miRNAs.
4. Prefer the sub-saturation graded calls (DIAPH1, ARID1B, CAMTA1) as the
   trustworthy examples of a real switch, independent of the miRNA story.

## 5. Gene track figures (geneviews) — named genes, by cell type

These are the PeakATail per-gene PAS tracks for named genes that were rendered
locally, one per (gene, cell type). Paths are relative to this file
(`reports/`). Note: only genes in the report's top-switch geneview selection were
pre-rendered — of the genes discussed above, **NUDT4, MSI2, AKAP13, NR3C1, FOXP1**
have local tracks; the rest (MDM2, HK1, GLS, CD44, CDKN1B, NF1, TIMP2, TGFBR2 …)
were **not** rendered and would need on-demand generation (the sweep produces
geneviews per run, or I can render specific ones from the cell-type `_inputs`
h5ads).

**NUDT4 — Macrophage M2 (shortening; subset_atlas)**
![NUDT4 Macrophage M2](laughney_full_report/figures/geneview/subset_atlas/CELL_TYPES_WANSLEEBEN_HOGAN_2013_MACROPHAGE_M2/figures/gene_ENSG00000173598.png)

**MSI2 — Macrophage M2 (lengthening, brain-met; subset_noatlas)**
![MSI2 Macrophage M2](laughney_full_report/figures/geneview/subset_noatlas/CELL_TYPES_WANSLEEBEN_HOGAN_2013_MACROPHAGE_M2/figures/gene_ENSG00000153944.png)

**AKAP13 — Mucinous epithelial (shortening; subset_atlas)**
![AKAP13 Mucinous](laughney_full_report/figures/geneview/subset_atlas/lung_epithelial_lineage_signatures_MUCINOUS/figures/gene_ENSG00000170776.png)

**NR3C1 / glucocorticoid receptor — AE2 (lengthening, validated miR-124/142-3p; b1_noatlas)**
![NR3C1 AE2](laughney_full_report/figures/geneview/b1_noatlas/CELL_TYPES_WANSLEEBEN_HOGAN_2013_AE2/figures/gene_ENSG00000113580.png)

**FOXP1 — Mucinous epithelial (shortening; b1_noatlas)**
![FOXP1 Mucinous](laughney_full_report/figures/geneview/b1_noatlas/lung_epithelial_lineage_signatures_MUCINOUS/figures/gene_ENSG00000114861.png)

Each has a `.svg` sibling (vector) and a `_pas_distances.csv` distance table in the
same folder.

---

## 6. References — internet sources used

**Mechanism / APA-in-cancer background**
- Mayr & Bartel 2009, *Cell* — 3′UTR shortening activates oncogenes: https://www.cell.com/fulltext/S0092-8674(09)00716-8 · https://pmc.ncbi.nlm.nih.gov/articles/PMC2819821/
- Park et al. 2018, *Nat Genet* — 3′UTR shortening represses tumor suppressors via ceRNA: https://www.nature.com/articles/s41588-018-0118-8
- LUAD APA landscape & TME: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8012674/
- Alternative polyadenylation shapes LUAD, *Hum Mol Genet* 2025 34(1):1: https://academic.oup.com/hmg/article/34/1/1/7866773

**Gene function / cancer role**
- MSI2 drives NSCLC metastasis, *PNAS* 2016: https://www.pnas.org/doi/10.1073/pnas.1513616113 · MSI2/EGFR *Oncogenesis* 2021: https://www.nature.com/articles/s41389-021-00317-y · MSI2 prognostic PMC8024834: https://pmc.ncbi.nlm.nih.gov/articles/PMC8024834/
- CD44 isoform switching & lung colonization, *Nat Commun* 2012: https://www.nature.com/articles/ncomms1892
- GLS/glutaminase targeting in NSCLC: https://pmc.ncbi.nlm.nih.gov/articles/PMC6150022/ · https://www.nature.com/articles/s12276-023-00971-9

**Validated miRNA–3′UTR interactions (the mechanism check)**
- CD44 / miR-34a (2 sites, luciferase): https://elifesciences.org/articles/06434 · renal https://pubmed.ncbi.nlm.nih.gov/31506763/ ; CD44 / miR-373-3p: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8230484/
- CDKN1B/p27 / miR-221–222 (EMBO J 2007): https://pubmed.ncbi.nlm.nih.gov/17627278/ · https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3567044/
- NF1 / miR-193b: https://pmc.ncbi.nlm.nih.gov/articles/PMC3546079/ · miR-107: https://www.nature.com/articles/srep36531 · miR-103: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3462785/
- TIMP2 / miR-197-3p (LUAD metastasis, luciferase): https://www.nature.com/articles/s41419-022-05420-5 · miR-210: https://pubmed.ncbi.nlm.nih.gov/27018975/
- TGFBR2 / miR-301a: https://link.springer.com/article/10.1186/s13046-014-0113-6 · miR-93: https://link.springer.com/article/10.1186/1476-4598-13-51 · miR-135b: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4462589/
- MECOM/EVI1 / miR-22: https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006259 · miR-133/Evi1: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4709720/
- FOXP1 / miR-34a: https://www.researchgate.net/publication/327042714_MicroRNA_miR-34a_downregulates_FOXP1
- NR3C1 / miR-124 (lengthening example): https://www.sciencedirect.com/science/article/abs/pii/S027858461730338X ; ATG7 / miR-17: https://pmc.ncbi.nlm.nih.gov/articles/PMC3742487/
- CERS6 3′UTR (predicted only) review: https://pmc.ncbi.nlm.nih.gov/articles/PMC6315813/

**Databases queried:** TargetScan Human (targetscan.org), miRTarBase (validated interactions), miRDB — used to check predicted vs experimentally-validated 3′UTR sites per gene. Gene oncogene/tumor-suppressor status reflects consensus (NCBI Gene / UniProt / literature above); APA-direction interpretations are inference from the Mayr/LUAD model, not gene-specific APA measurements in the cited papers.
