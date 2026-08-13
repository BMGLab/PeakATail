#!/usr/bin/env bash
# Rebuild the LaTeX technical report (two passes for TOC + cross-references).
set -euo pipefail
# Resolve the report directory BEFORE any cd so $(dirname ...) doesn't drift.
REPORTS_DIR="$(cd "$(dirname "$0")" && pwd)"
# A fresh directory per build. Reusing one carries .aux/.toc across runs, so
# pdflatex converges from the previous run's page numbering instead of from
# scratch -- that made the reported page count drift between edits (48 vs 45 for
# identical content) and masked a section having moved into the appendix.
BUILDDIR=$(mktemp -d "${TMPDIR:-/tmp}/latex_build.XXXXXX")
trap 'rm -f "$BUILDDIR"/*.aux "$BUILDDIR"/*.toc "$BUILDDIR"/*.out 2>/dev/null || true' EXIT
cp -r "$REPORTS_DIR/figures" "$BUILDDIR/"
cp "$REPORTS_DIR/PeakATail_Technical_Report.tex" "$BUILDDIR/"
cd "$BUILDDIR"
pdflatex -interaction=nonstopmode PeakATail_Technical_Report.tex > build1.log 2>&1
pdflatex -interaction=nonstopmode PeakATail_Technical_Report.tex > build2.log 2>&1

# Fail loudly: a nonstopmode build "succeeds" even when figures are missing or
# references dangle, so the report would ship with blank boxes unnoticed.
if grep -qE '^!|Unable to load' build2.log; then
    echo "ERROR: LaTeX reported errors:" >&2
    grep -nE '^!|Unable to load' build2.log | head -20 >&2
    exit 1
fi
if grep -q 'undefined' build2.log; then
    echo "WARNING: undefined references:" >&2
    grep -n 'undefined on input line' build2.log | head -10 >&2
fi
cp PeakATail_Technical_Report.pdf "$REPORTS_DIR/"
pages=$(pdfinfo PeakATail_Technical_Report.pdf 2>/dev/null | awk '/^Pages/{print $2}')
figs=$(tr -d '\n' < build2.log | grep -o 'figures/corrected/[0-9]*_[a-z_]*\.png' | sort -u | wc -l)
echo "Build complete: $REPORTS_DIR/PeakATail_Technical_Report.pdf (${pages} pages, ${figs} corrected figures)"
