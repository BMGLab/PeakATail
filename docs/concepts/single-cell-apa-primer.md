# Single-cell APA primer

This page explains alternative polyadenylation (APA) for a reader who is comfortable with single-cell RNA-seq but has not studied 3' end biology. It is background reading before you interpret PeakATail results.

---

## What is alternative polyadenylation?

Every protein-coding messenger RNA (mRNA) ends with a poly(A) tail — a stretch of adenosine nucleotides added after the stop codon. The site where the tail is added is called a **poly(A) site** (PAS). Most genes in the human genome have more than one PAS; the cell can choose which one to use when it processes the pre-mRNA. When different cells or conditions use different PAS within the same gene, the result is **alternative polyadenylation**.

The critical consequence is 3'UTR length. The 3' untranslated region (3'UTR) is the segment of mRNA between the stop codon and the poly(A) tail. A PAS that lies close to the stop codon (a *proximal* PAS) produces a short 3'UTR; a PAS far from the stop codon (a *distal* PAS) produces a long 3'UTR. The coding sequence is identical in both cases, but the regulatory content of the transcript is not.

Long 3'UTRs contain binding sites for microRNAs, RNA-binding proteins, and other post-transcriptional regulators. Shifting from a distal to a proximal PAS removes these binding sites, potentially increasing protein output, altering mRNA stability, or changing subcellular localisation. APA is therefore a mechanism for post-transcriptional regulation that is independent of transcriptional changes.

---

## Why does APA matter in single-cell biology?

Bulk RNA-seq averages APA signals over thousands of cells and can detect only the most extreme shifts between conditions. Single-cell 3' end sequencing reveals that APA usage is **cell-type specific**: a T cell and a monocyte from the same tissue sample may express the same gene at the same level but use entirely different PAS. This cell-type specificity is invisible to bulk profiling and is missed by most single-cell pipelines that focus on gene-level counts.

Several published studies have documented this:

- Bai et al. (2020) showed cell-type-specific APA in complex tissues using single-cell data, finding that APA patterns distinguish cell types that are otherwise similar by gene expression alone. (DOI: [10.1101/2020.07.30.229096](https://doi.org/10.1101/2020.07.30.229096))
- Tian et al. (2005) established the role of APA in regulating 3'UTR length across tissues and development stages, providing the foundational evidence that proximal/distal PAS selection is a widespread regulatory layer. (DOI: [10.1038/ng1533](https://doi.org/10.1038/ng1533))
- Wu et al. (2021) developed scAPAtrap, demonstrating that single-cell APA sites can be detected genome-wide and that APA usage heterogeneity between cell types is a reproducible biological signal distinct from gene expression differences. (DOI: [10.1093/bib/bbaa273](https://doi.org/10.1093/bib/bbaa273))

---

## How does 3' end sequencing detect PAS?

10x Chromium scRNA-seq is inherently a 3' end protocol: reverse transcription is primed from the poly(A) tail, so sequenced reads preferentially capture the 3' end of each transcript. This means that read pile-ups in the BAM file represent genuine polyadenylation sites, not random sampling of the transcript body.

PeakATail exploits this. It scans the BAM file strand by strand, counts how many read 3' ends fall at each genomic position, and applies a peak-calling algorithm to identify positions with statistically elevated coverage (the actual PAS). The result is a set of PAS coordinates — one BED entry per detected site — and a count matrix recording how many reads from each cell barcode ended at each PAS.

The matrix is sparse: most cells have zero reads at most PAS. This is expected and analogous to the sparsity of gene expression matrices in scRNA-seq. PeakATail's preprocessing applies the same filtering logic (minimum reads per cell, minimum cells per PAS) that Seurat and Scanpy apply to gene expression matrices.

---

## What do proximal, distal, and PDUI mean?

Given a gene with two or more detected PAS:

- The **proximal PAS** is the one closest to the stop codon (smallest genomic coordinate on the + strand; largest coordinate on the − strand). It produces the shortest 3'UTR.
- The **distal PAS** is the one farthest from the stop codon. It produces the longest 3'UTR.

For a gene with exactly two PAS, the **Proximal-Distal Usage Index** (PDUI) is:

```
PDUI = reads_at_distal / (reads_at_proximal + reads_at_distal)
```

PDUI ranges from 0 (all reads at the proximal PAS, short 3'UTR) to 1 (all reads at the distal PAS, long 3'UTR). A PDUI of 0.5 means equal usage of both isoforms.

A **delta-PDUI** between two cell clusters is the difference in their mean PDUI:

- Positive delta-PDUI: cells in the comparison group use longer 3'UTR isoforms (3'UTR lengthening).
- Negative delta-PDUI: cells in the reference group use longer 3'UTR isoforms (3'UTR shortening in the comparison group).

In activated immune cells, for example, 3'UTR shortening (shift toward proximal PAS) is frequently observed and is thought to increase protein output by removing microRNA-mediated suppression. PeakATail's `peakatail switch length --strategy classic` computes this PDUI value per cell per gene.

For genes with more than two PAS, the `proportion` strategy reports the fraction of reads at each individual PAS (not just the proximal/distal pair), and the `shannon` strategy reports the Shannon entropy of the per-PAS distribution — higher entropy = more uniform usage across multiple PAS.
