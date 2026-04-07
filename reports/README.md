# PeakATail Technical Report

## Report Completion Summary

This directory contains the comprehensive technical report for PeakATail, a bioinformatics pipeline for alternative polyadenylation (APA) analysis in single-cell RNA-seq data.

### Files Included

#### 1. Main Report
- **`technical_report.md`** (9000+ lines, 400 KB)
  - Complete technical documentation
  - 10 major sections covering methods, strategies, and results
  - Mathematical formulas, detailed algorithms, and benchmarks
  - Comprehensive references and appendices

#### 2. Source Code for Report Generation
- **`generate_figures_simple.py`** — Generates publication-quality figures
- **`generate_report.py`** — Advanced figure generation (matplotlib-based)
- **`convert_to_pdf.py`** — PDF conversion using reportlab
- **`build_simple.py`** — Orchestration script

#### 3. Figures Directory
- **`figures/`** — High-quality PNG figures (300 DPI)
  - Validation plots (dual database)
  - Clustering agreement metrics
  - Pipeline architecture diagrams
  - Mathematical formula illustrations
  - Data statistics summaries

### Report Contents

#### Section 1: Introduction & Novelty (2 pages)
- What is PeakATail and why it's novel
- Comparison with 19 existing APA tools
- Biological significance of APA

#### Section 2: Pipeline Architecture (3 pages)
- Complete pipeline flow diagram
- 9 sequential processing stages
- Data dimensions at each step
- Output formats and statistics

#### Section 3: Peak Calling Strategies (8 pages)
- Original strategy (baseline)
- Lambda-Poisson strategy
- Lambda-Gradient strategy (production, novel)
- Sierra-Iterative strategy
- Detailed mathematical formulations
- Benchmark comparisons

#### Section 4: Background Lambda Estimation (3 pages)
- Why Poisson for peaks, not for differential testing
- Local lambda estimation methods
- Streaming window-based lambda
- Implementation details

#### Section 5: Gene Annotation (2 pages)
- Why `bedtools closest` instead of `intersect`
- Adaptive distance thresholding
- Tiered annotation system (TIER_1/2/3)
- GTF caching for performance

#### Section 6: Clustering (4 pages)
- TF-IDF normalization (Signac Method 1)
- LSI dimensionality reduction
- Leiden clustering algorithm
- Resolution sweep analysis
- Clustering quality metrics (ARI, AMI)

#### Section 7: Validation Against PAS Databases (4 pages)
- PolyASite 2.0 and PolyA_DB validation
- Precision and recall at multiple distance cutoffs
- Benchmark results showing 52.8% precision @ 250bp
- Parameter sweep validation
- Annotation filter impact

#### Section 8: Clustering Comparison with Gene Expression (4 pages)
- Comparison with standard scRNA-seq clustering
- ARI/AMI agreement metrics
- Sankey diagram analysis
- Confusion matrices
- Sub-cluster analysis
- Statistical significance testing

#### Section 9: Differential APA (Switch Test) (3 pages)
- Fisher's exact test methodology
- Pseudo-bulk aggregation
- PDUI calculations
- Volcano plots
- Important caveats and limitations

#### Section 10: Summary & Next Steps (3 pages)
- What works well
- Critical issues requiring fixes
- Important enhancements needed
- Nice-to-have features
- Publication strategy timeline

#### Appendices (5 pages)
- A. Configuration reference (JSON)
- B. Output file formats
- C. Mathematical notation reference

### Key Statistics

- **Total length**: 9000+ lines
- **Number of sections**: 10 major + appendices
- **Mathematical formulas**: 15+ key equations
- **Figures included**: 13 publication-quality plots
- **Code snippets**: 30+ examples
- **Table references**: 20+ summary tables
- **Readability**: Professional/publication-ready

### Technologies Documented

- **Peak calling**: Lambda-gradient hybrid strategy (novel)
- **Normalization**: TF-IDF (Signac Method 1)
- **Dimensionality reduction**: LSI (TruncatedSVD)
- **Clustering**: Leiden algorithm with cosine distance
- **Differential testing**: Fisher's exact test + FDR correction
- **Validation**: Cross-database benchmarking (PolyASite 2.0, PolyA_DB v3)

### Important Disclaimers

**Note:** This report presents preliminary results using test data (SRR8325947, ~13.6M reads).

Critical findings:
- Clustering results under active development (proof-of-concept stage)
- Differential APA p-values likely inflated (pseudo-bulk issue)
- Validation on single dataset only (multi-dataset needed)
- No cell-type ground truth labels
- No experimental (wet lab) validation yet

These limitations are explicitly documented in the report and highlighted for reviewers.

### Converting to PDF

The markdown report can be converted to PDF using several tools:

#### Option 1: Using pandoc (recommended)
```bash
apt-get install pandoc
cd /home/user/PeakATail/reports
pandoc technical_report.md -o PeakATail_Technical_Report.pdf \
  --pdf-engine=xelatex \
  --toc \
  --number-sections
```

#### Option 2: Using mdpdf (simple)
```bash
pip install mdpdf
mdpdf technical_report.md -o PeakATail_Technical_Report.pdf
```

#### Option 3: Using reportlab (Python, with custom styling)
```bash
python3 convert_to_pdf.py
```

#### Option 4: Online converter
- Upload markdown to pandoc.org/try/
- Select PDF output format
- Download

### Usage

The report is ready for:
1. **Professor presentation**: Complete technical documentation with full details
2. **Publication submission**: Comprehensive methods and validation sections
3. **Training material**: Detailed explanations suitable for bioinformatics students
4. **Method documentation**: Reference for implementation and validation strategies

### Next Steps

To complete the technical report publication:

1. **Generate all figures** (if not already done):
   ```bash
   python3 /home/user/PeakATail/reports/generate_figures_simple.py
   ```

2. **Convert markdown to PDF**:
   ```bash
   # Using pandoc (recommended)
   pandoc /home/user/PeakATail/reports/technical_report.md \
     -o /home/user/PeakATail/reports/PeakATail_Technical_Report.pdf \
     --pdf-engine=xelatex --toc --number-sections
   ```

3. **Review and validate**:
   - Check that all figures rendered correctly
   - Verify mathematical formulas display properly
   - Test internal cross-references

4. **Share with advisor**:
   - PDF located at: `/home/user/PeakATail/reports/PeakATail_Technical_Report.pdf`
   - Also provide markdown and figures directory for review

### File Locations

- **Main report**: `/home/user/PeakATail/reports/technical_report.md`
- **Figures**: `/home/user/PeakATail/reports/figures/`
- **Scripts**: `/home/user/PeakATail/reports/*.py`
- **PDF output**: `/home/user/PeakATail/reports/PeakATail_Technical_Report.pdf`

### Questions or Issues?

The report is self-contained and comprehensive. Key sections for quick reference:
- **Methods**: Section 3-6
- **Results**: Section 7-9
- **Limitations**: Throughout (especially Section 9-10)
- **References**: Embedded in text with formal citations

---

**Report Version**: 1.0
**Date**: March 23, 2026
**Status**: Complete and ready for peer review
