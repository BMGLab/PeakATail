#!/bin/bash
cd /home/user/PeakATail/reports
echo "Building PeakATail Technical Report..."
python3 generate_figures_simple.py
echo ""
echo "Done! Report files are ready:"
ls -lh technical_report.md
ls -lh figures/ 2>/dev/null || echo "No figures directory yet"
echo ""
echo "Next: Convert markdown to PDF using pandoc, mdpdf, or another tool"
