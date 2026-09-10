# Concepts

These pages explain what is happening inside PeakATail before you run it. Read them when you want to understand *why* the pipeline is structured the way it is, not just *how* to invoke it.

## Before you run — conceptual map

| Question | Where to look |
|----------|--------------|
| How does data move from BAM to clusters.h5ad? | [Data flow](data-flow.md) |
| What is every file in my output directory? | [Output files](output-files.md) |
| What is APA and why does it matter in single-cell data? | [Single-cell APA primer](single-cell-apa-primer.md) |

## Page descriptions

**[Data flow](data-flow.md)**
A Mermaid diagram of every pipeline stage with the file types that flow between them. Covers both the single-sample path (one BAM) and the multi-sample path (two or more BAMs with atlas snap or coordinate merge).

**[Output files](output-files.md)**
A reference entry for every file PeakATail writes. Includes the canonical path, format, and the pipeline stage that produces it. Useful when you want to load a specific artifact into Python or pass it to another tool.

**[Single-cell APA primer](single-cell-apa-primer.md)**
A ~600-word explanation of alternative polyadenylation for a reader who knows single-cell sequencing but has not worked with 3' end biology before. Covers what a PAS is, why 3'UTR length is regulated per cell type, and how PDUI quantifies it.
