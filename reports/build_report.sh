#!/bin/bash
set -e

echo "=========================================================================="
echo "PeakATail Technical Report: Complete Build"
echo "=========================================================================="

cd /home/user/PeakATail/reports

# Step 1: Generate figures
echo ""
echo "Step 1: Generating matplotlib figures..."
python3 generate_report.py

# Step 2: Convert to PDF
echo ""
echo "Step 2: Converting markdown + figures to PDF..."
python3 convert_to_pdf.py

# Step 3: Summary
echo ""
echo "=========================================================================="
echo "Build complete!"
echo "=========================================================================="
ls -lh PeakATail_Technical_Report.pdf
echo ""
echo "Report saved to: /home/user/PeakATail/reports/PeakATail_Technical_Report.pdf"
echo "Figures saved to: /home/user/PeakATail/reports/figures/"
