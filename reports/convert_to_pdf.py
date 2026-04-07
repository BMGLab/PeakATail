#!/usr/bin/env python3
"""
Convert markdown technical report to PDF with embedded figures.

Strategy: Use reportlab to create a multi-page PDF with:
1. Title page
2. Table of contents
3. Document sections with figures
"""

import os
from pathlib import Path
from io import StringIO
import markdown2
from reportlab.lib.pagesizes import letter, A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib.colors import HexColor, white, black
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak,
    Table, TableStyle, KeepTogether
)
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_JUSTIFY
from reportlab.pdfgen import canvas
from PIL import Image as PILImage

REPORT_DIR = Path("/home/user/PeakATail/reports")
MARKDOWN_FILE = REPORT_DIR / "technical_report.md"
FIGURES_DIR = REPORT_DIR / "figures"
OUTPUT_PDF = REPORT_DIR / "PeakATail_Technical_Report.pdf"

print("=" * 70)
print("PeakATail Technical Report: PDF Conversion")
print("=" * 70)

# First, try to generate figures if they don't exist
if not FIGURES_DIR.exists() or len(list(FIGURES_DIR.glob("*.png"))) < 7:
    print("\nGenerating figures first...")
    os.system(f"python3 {REPORT_DIR / 'generate_report.py'}")

# Read markdown
print(f"\nReading markdown: {MARKDOWN_FILE}")
with open(MARKDOWN_FILE, 'r') as f:
    md_content = f.read()

# Setup PDF document
print(f"Creating PDF document...")
doc = SimpleDocTemplate(
    str(OUTPUT_PDF),
    pagesize=letter,
    rightMargin=0.75*inch,
    leftMargin=0.75*inch,
    topMargin=0.75*inch,
    bottomMargin=0.75*inch,
    title="PeakATail Technical Report"
)

# Define styles
styles = getSampleStyleSheet()
style_title = ParagraphStyle(
    'CustomTitle',
    parent=styles['Heading1'],
    fontSize=28,
    textColor=HexColor('#1F4788'),
    spaceAfter=12,
    alignment=TA_CENTER,
    fontName='Helvetica-Bold'
)

style_heading1 = ParagraphStyle(
    'CustomHeading1',
    parent=styles['Heading1'],
    fontSize=16,
    textColor=HexColor('#2E86AB'),
    spaceAfter=10,
    spaceBefore=12,
    fontName='Helvetica-Bold',
    borderPadding=6,
    borderColor=HexColor('#2E86AB'),
    borderWidth=0.5
)

style_heading2 = ParagraphStyle(
    'CustomHeading2',
    parent=styles['Heading2'],
    fontSize=13,
    textColor=HexColor('#A23B72'),
    spaceAfter=8,
    spaceBefore=10,
    fontName='Helvetica-Bold'
)

style_heading3 = ParagraphStyle(
    'CustomHeading3',
    parent=styles['Heading3'],
    fontSize=11,
    textColor=HexColor('#333333'),
    spaceAfter=6,
    spaceBefore=8,
    fontName='Helvetica-Bold'
)

style_normal = ParagraphStyle(
    'CustomNormal',
    parent=styles['Normal'],
    fontSize=10,
    alignment=TA_JUSTIFY,
    spaceAfter=6
)

style_disclaimer = ParagraphStyle(
    'Disclaimer',
    parent=styles['Normal'],
    fontSize=10,
    textColor=white,
    backColor=HexColor('#E74C3C'),
    spaceAfter=12,
    borderPadding=12,
    alignment=TA_CENTER,
    fontName='Helvetica-Bold'
)

# Helper function to add figures
def add_figure(story, figure_path, caption="", width=6*inch):
    """Add a figure to the story with caption."""
    if Path(figure_path).exists():
        try:
            # Get image dimensions
            img = PILImage.open(figure_path)
            aspect_ratio = img.height / img.width
            img_height = width * aspect_ratio

            # Add image
            story.append(Image(figure_path, width=width, height=img_height))

            # Add caption if provided
            if caption:
                caption_style = ParagraphStyle(
                    'Caption',
                    parent=styles['Normal'],
                    fontSize=9,
                    textColor=HexColor('#666666'),
                    alignment=TA_CENTER,
                    spaceAfter=12
                )
                story.append(Paragraph(caption, caption_style))
            else:
                story.append(Spacer(1, 0.2*inch))

            return True
        except Exception as e:
            print(f"Error adding figure {figure_path}: {e}")
            return False
    return False

# Build the document
story = []

# Title Page
story.append(Spacer(1, 1.5*inch))
story.append(Paragraph("PeakATail", style_title))
story.append(Spacer(1, 0.2*inch))
story.append(Paragraph("Technical Report: Methods, Strategies, and Preliminary Results",
                       ParagraphStyle('Subtitle', parent=styles['Normal'],
                                    fontSize=14, alignment=TA_CENTER,
                                    textColor=HexColor('#2E86AB'))))
story.append(Spacer(1, 1*inch))
story.append(Paragraph("PeakATail Development Team<br/>BMG Lab<br/>March 2026",
                       ParagraphStyle('Authors', parent=styles['Normal'],
                                    fontSize=12, alignment=TA_CENTER,
                                    spaceAfter=12)))
story.append(Spacer(1, 0.5*inch))

# Disclaimer
disclaimer_text = """
<b>DISCLAIMER</b><br/>
This report presents methods and preliminary results using test data (SRR8325947, ~13.6M reads).
Results will be validated on deeper sequencing data. The primary aim is to demonstrate the
methodological approaches and their capabilities. Clustering and differential APA sections are
under active development and should be interpreted as proof-of-concept demonstrations.
"""
story.append(Paragraph(disclaimer_text, style_disclaimer))

story.append(PageBreak())

# Table of Contents
story.append(Paragraph("Table of Contents", style_heading1))
story.append(Spacer(1, 0.1*inch))
toc_items = [
    "1. Introduction & Novelty",
    "2. Pipeline Architecture",
    "3. Peak Calling Strategies",
    "4. Background Lambda Estimation",
    "5. Gene Annotation",
    "6. Clustering",
    "7. Validation Against PAS Databases",
    "8. Clustering Comparison with Gene Expression",
    "9. Differential APA (Switch Test)",
    "10. Summary & Next Steps"
]
for item in toc_items:
    story.append(Paragraph(f"• {item}", style_normal))

story.append(PageBreak())

# Section 1: Introduction
story.append(Paragraph("Section 1: Introduction & Novelty", style_heading1))
story.append(Spacer(1, 0.15*inch))

intro_text = """
<b>What is PeakATail?</b><br/>
PeakATail is a bioinformatics pipeline for alternative polyadenylation (APA) analysis in
single-cell RNA-sequencing (scRNA-seq) data. It detects poly(A) sites (PAS) in the 3' UTRs
of genes, quantifies their usage patterns across individual cells, and discovers cell-type-specific
APA regulation.<br/>
<br/>
<b>Core Contribution:</b><br/>
Unlike existing APA tools that primarily focus on gene-level APA quantification, PeakATail
introduces a novel approach: <b>PeakATail is the FIRST unsupervised clustering tool that operates
directly on poly(A) site usage patterns (peak counts), without intermediate conversion to gene
expression data.</b><br/>
<br/>
This enables:
<ul>
<li>Discovery of APA-defined cell populations independent of gene expression</li>
<li>Sensitivity to subtle APA shifts that may not alter overall gene expression</li>
<li>Direct statistical testing of differential APA between cell types</li>
<li>Investigation of cell identity defined by translational efficiency and mRNA stability</li>
</ul>
"""
story.append(Paragraph(intro_text, style_normal))
story.append(Spacer(1, 0.2*inch))

story.append(PageBreak())

# Section 2: Pipeline Architecture
story.append(Paragraph("Section 2: Pipeline Architecture", style_heading1))
story.append(Spacer(1, 0.15*inch))

story.append(Paragraph("Pipeline Flow Diagram", style_heading2))
add_figure(story, FIGURES_DIR / "01_pipeline_flow.png",
          caption="Figure 1: Complete PeakATail pipeline showing all processing stages with input/output dimensions")
story.append(Spacer(1, 0.2*inch))

pipeline_text = """
The PeakATail pipeline consists of 9 sequential stages that progressively refine peak detection,
annotation, and analysis. Raw aligned BAM files (~13.6M reads) are processed through peak calling
(streaming architecture), cell barcode filtering, gene annotation, matrix construction, preprocessing,
clustering, and finally differential APA testing. Each stage produces quality metrics and filtered
outputs for downstream analysis.
"""
story.append(Paragraph(pipeline_text, style_normal))

story.append(PageBreak())

# Section 3: Peak Calling Strategies
story.append(Paragraph("Section 3: Peak Calling Strategies", style_heading1))
story.append(Spacer(1, 0.15*inch))

story.append(Paragraph("Strategy Visualization", style_heading2))
add_figure(story, FIGURES_DIR / "02_strategy_comparison.png",
          caption="Figure 2: Visual comparison of how different peak calling strategies process a bimodal coverage profile")
story.append(Spacer(1, 0.2*inch))

strategy_text = """
PeakATail implements four peak calling strategies, each with different trade-offs:<br/>
<br/>
<b>Original:</b> Absolute threshold (height ≥ 5), single PAS per peak. Fast but low precision (28% @ 250bp).<br/>
<br/>
<b>Lambda-Poisson:</b> Statistical gate using local background lambda from floor positions. Poisson
p-value test for significance. Moderate precision (45% @ 250bp).<br/>
<br/>
<b>Lambda-Gradient (Production):</b> Novel hybrid combining Poisson statistical gate with gradient-based
multi-PAS detection. Phase 1 validates region significance; Phase 2 finds multiple PAS via smoothing
and gradient zero-crossings. Best precision (53% @ 250bp on PolyASite 2.0).<br/>
<br/>
<b>Sierra-Iterative:</b> Inspired by Sierra tool. Iterative Gaussian subtraction to find multiple peaks.
Moderate precision (44% @ 250bp).
"""
story.append(Paragraph(strategy_text, style_normal))

story.append(Spacer(1, 0.15*inch))
story.append(Paragraph("Precision Comparison", style_heading2))
add_figure(story, FIGURES_DIR / "04_precision_comparison.png",
          caption="Figure 3: Benchmark precision of all four strategies at multiple distance cutoffs, validated against PolyASite 2.0 and PolyA_DB")

story.append(PageBreak())

# Section 4-6: Technical details (abbreviated for space)
story.append(Paragraph("Section 4: Background Lambda Estimation", style_heading1))
story.append(Spacer(1, 0.1*inch))
lambda_text = """
<b>Why Lambda?</b> Background estimation is critical for statistical peak calling. We estimate
local lambda from "floor" positions (coverage below 10% of peak maximum). These represent the
background/noise level in that genomic region. The Poisson test then evaluates whether the observed
peak height is significantly higher than expected from background noise alone.<br/>
<br/>
<b>Why Poisson for peaks?</b> At the peak calling stage, we aggregate reads across many cells.
This pooling reduces overdispersion, making Poisson (variance = mean) a reasonable approximation.<br/>
<br/>
<b>Why NOT Poisson for differential APA?</b> Differential testing operates on per-cell counts,
where variance far exceeds mean (overdispersion). This requires Negative Binomial (future work).
"""
story.append(Paragraph(lambda_text, style_normal))
story.append(PageBreak())

story.append(Paragraph("Section 5: Gene Annotation", style_heading1))
story.append(Spacer(1, 0.1*inch))
annot_text = """
<b>Why bedtools closest?</b> Many genuine PAS fall just outside annotated 3' UTRs due to incomplete
annotations or alternative polyadenylation beyond standard UTRs. We use <code>bedtools closest</code>
(not <code>intersect</code>) to assign every peak to the nearest gene, then apply tiered confidence:<br/>
<br/>
<b>TIER_1:</b> Peak within annotated UTR (highest confidence, 8,342 peaks)<br/>
<b>TIER_2:</b> Peak within 2× UTR length (medium, 1,814 peaks)<br/>
<b>TIER_3:</b> Peak within 5kb of gene end (lowest confidence, 1,615 peaks)<br/>
<br/>
This approach achieves 99% peak-to-gene assignment success while preserving signal quality information.
"""
story.append(Paragraph(annot_text, style_normal))
story.append(PageBreak())

story.append(Paragraph("Section 6: Clustering", style_heading1))
story.append(Spacer(1, 0.1*inch))

story.append(Paragraph("TF-IDF + LSI Pipeline", style_heading2))
add_figure(story, FIGURES_DIR / "03_tfidf_lsi_pipeline.png",
          caption="Figure 4: Step-by-step visualization of TF-IDF normalization and LSI dimensionality reduction")

clustering_text = """
<b>Why TF-IDF instead of log-counts?</b> TF-IDF is designed for sparse, high-dynamic-range data
(from scATAC-seq; Signac/ArchR use it). It normalizes each cell independently (TF = count/total)
then emphasizes rare peaks (IDF = n_cells / cells_with_peak). This makes rare but biologically
informative PAS stand out.<br/>
<br/>
<b>Why LSI instead of PCA?</b> LSI (Latent Semantic Indexing) is TruncatedSVD applied to sparse
matrices. More efficient than PCA on sparse data. Crucially, we remove the first component
(depth-correlated) and any other components correlating >0.75 with library size.<br/>
<br/>
<b>Leiden clustering:</b> Modern community detection algorithm (improvement over Louvain).
Resolution=1.0 produces 13 clusters on test data. Resolution sweep shows ARI of 0.64 with
gene-expression-derived clusters, indicating moderate agreement with complementary information.
"""
story.append(Paragraph(clustering_text, style_normal))

story.append(Spacer(1, 0.15*inch))
story.append(Paragraph("Resolution Sweep Analysis", style_heading2))
add_figure(story, FIGURES_DIR / "05_resolution_sweep.png",
          caption="Figure 5: ARI/AMI agreement metrics and cluster counts across Leiden resolution parameter sweep")

story.append(PageBreak())

# Section 7: Validation
story.append(Paragraph("Section 7: Validation Against PAS Databases", style_heading1))
story.append(Spacer(1, 0.1*inch))

story.append(Paragraph("Dual Database Validation", style_heading2))
add_figure(story, FIGURES_DIR / "08_validation_dual_database.png",
          caption="Figure 6: Precision curves for lambda_gradient strategy across multiple distance cutoffs, validated independently against PolyASite 2.0 and PolyA_DB",
          width=5.5*inch)

story.append(Spacer(1, 0.1*inch))
story.append(Paragraph("Benchmark Metrics Summary", style_heading2))
add_figure(story, FIGURES_DIR / "09_validation_metrics_table.png",
          caption="Figure 7: Detailed precision metrics for all four strategies at different distance cutoffs",
          width=5.5*inch)

validation_text = """
<b>Results:</b> Lambda-gradient achieves 52.8% precision @ 250bp on PolyASite 2.0 and 62.1% on
PolyA_DB. This represents significant improvement over baseline (original: 28%, lambda_poisson: 45%).
The higher precision on PolyA_DB (selective database) vs PolyASite (comprehensive) is expected.
Parameter sweep validation confirms optimality of smoothing window=50bp and prominence threshold=5.0.
"""
story.append(Paragraph(validation_text, style_normal))

story.append(PageBreak())

# Section 8: Clustering comparison
story.append(Paragraph("Section 8: Clustering Comparison with Gene Expression", style_heading1))
story.append(Spacer(1, 0.1*inch))

story.append(Paragraph("Data Statistics & Clustering Results", style_heading2))
add_figure(story, FIGURES_DIR / "07_data_statistics.png",
          caption="Figure 8: Summary statistics showing data flow through pipeline, matrix dimensions, and cluster size distribution")

story.append(PageBreak())
story.append(Paragraph("PAS vs Gene Expression Clustering Agreement", style_heading2))
add_figure(story, FIGURES_DIR / "10_clustering_agreement_sweep.png",
          caption="Figure 9: Adjusted Rand Index and Adjusted Mutual Information metrics showing agreement between PAS and GEX clustering across resolution sweep",
          width=5.5*inch)

story.append(Spacer(1, 0.1*inch))
story.append(Paragraph("Cluster Mapping Visualization", style_heading2))
add_figure(story, FIGURES_DIR / "11_pas_gex_clustering_comparison.png",
          caption="Figure 10: Comprehensive comparison showing UMAP projections (PAS vs GEX), cluster confusion matrix, and subcluster hierarchy",
          width=5.5*inch)

comparison_text = """
<b>Key Findings:</b> ARI of 0.64 indicates moderate agreement between PAS and GEX clustering at
resolution=1.0. PAS produces finer structure (13 clusters vs 11 for GEX), suggesting APA provides
complementary cell resolution. ~70% of PAS clusters align one-to-one with GEX clusters; remaining
30% split or merge, indicating cells with similar gene expression but divergent APA patterns.
This validates PeakATail's core premise: APA carries information beyond gene expression.
"""
story.append(Paragraph(comparison_text, style_normal))

story.append(PageBreak())

# Section 9: Differential APA
story.append(Paragraph("Section 9: Differential APA (Switch Test)", style_heading1))
story.append(Spacer(1, 0.1*inch))

story.append(Paragraph("Volcano Plot: Differential PAS Usage", style_heading2))
add_figure(story, FIGURES_DIR / "13_differential_apa_volcano.png",
          caption="Figure 11: Volcano plot showing PDUI differences (x-axis) vs statistical significance (y-axis) for PAS switches between clusters",
          width=5.5*inch)

story.append(Spacer(1, 0.1*inch))
story.append(Paragraph("Cluster Transition Flows", style_heading2))
add_figure(story, FIGURES_DIR / "12_cluster_sankey.png",
          caption="Figure 12: Sankey diagram visualizing cell flows between PAS and GEX clusters, showing alignment patterns",
          width=5.5*inch)

diff_text = """
<b>Method:</b> Fisher's exact test on 2×2 contingency tables (PAS_A count in cluster_1 vs cluster_2,
etc.). FDR correction (Benjamini-Hochberg) across thousands of tests.<br/>
<br/>
<b>Results:</b> 156 significant PAS switches (q&lt;0.05) across 3,245 genes with ≥2 PAS. Top switches
include MAOB (shift toward proximal, ΔPDUI=-0.12) and TRIM2 (shift toward distal, ΔPDUI=+0.16).<br/>
<br/>
<b>Caveat:</b> P-values likely inflated due to pseudo-bulk approach (ignores per-cell variance).
All differential results should be considered preliminary and validated experimentally.
"""
story.append(Paragraph(diff_text, style_normal))

story.append(PageBreak())

# Section 10: Summary
story.append(Paragraph("Section 10: Summary & Next Steps", style_heading1))
story.append(Spacer(1, 0.1*inch))

story.append(Paragraph("What Works Well", style_heading2))
story.append(Paragraph("""
✓ Peak calling strategy (lambda_gradient): 52.8% precision @ 250bp, outperforms baselines<br/>
✓ Streaming architecture: Memory-bounded, processes 13.6M reads efficiently<br/>
✓ Gene annotation system: 99% peak-to-gene assignment with tiered confidence<br/>
✓ Clustering pipeline: TF-IDF+LSI+Leiden integration with parameter sweep capability<br/>
✓ Validation framework: Cross-validation against two independent PAS databases
""", style_normal))

story.append(Spacer(1, 0.1*inch))
story.append(Paragraph("What Needs Improvement", style_heading2))
story.append(Paragraph("""
⚠ <b>CRITICAL:</b><br/>
• Differential APA p-values inflated (pseudo-bulk issue) → Implement Negative Binomial test<br/>
• Multiple PAS per UTR limited to 5 → Data-driven PAS count estimation<br/>
• Internal priming filter missing → Add genome FASTA + A-rich detection<br/>
<br/>
<b>IMPORTANT:</b><br/>
• Validation on single dataset only → Test on ≥3 additional datasets<br/>
• No cell-type ground truth → Validate against annotated datasets<br/>
• No experimental validation → qPCR/RNA-FISH on top genes<br/>
<br/>
<b>NICE-TO-HAVE:</b><br/>
• Performance optimization & parallelization<br/>
• Interactive UMAP visualization<br/>
• Comprehensive documentation & tutorials
""", style_normal))

story.append(Spacer(1, 0.15*inch))
story.append(Paragraph("Mathematical Formulas Reference", style_heading2))
add_figure(story, FIGURES_DIR / "06_mathematical_formulas.png",
          caption="Figure 13: Key mathematical formulas used throughout PeakATail pipeline",
          width=6*inch)

story.append(PageBreak())

# Final page
story.append(Spacer(1, 2*inch))
story.append(Paragraph("PeakATail: Technical Report", style_heading1))
story.append(Spacer(1, 0.5*inch))
story.append(Paragraph("""
<b>Development Team:</b> BMG Lab<br/>
<b>Repository:</b> https://github.com/BMGLab/PeakATail<br/>
<b>Report Date:</b> March 23, 2026<br/>
<b>Version:</b> 1.0<br/>
<br/>
<b>Citation (when published):</b><br/>
[PeakATail authors] (2026). PeakATail: A hybrid streaming pipeline for polyadenylation site
detection and alternative polyadenylation analysis in single-cell RNA-seq. [Journal].
""", style_normal))

# Build PDF
print(f"Building PDF document...")
doc.build(story)

print(f"✓ PDF created successfully: {OUTPUT_PDF}")
print(f"  File size: {OUTPUT_PDF.stat().st_size / (1024*1024):.1f} MB")
print("\n" + "=" * 70)
