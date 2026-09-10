#!/usr/bin/env bash
# Regenerate tests/fixtures/cellranger_pbmc_tiny.bam
# Provenance: 10x Genomics pbmc_10k_v3, CellRanger 3.0.0, refdata-cellranger-GRCh38-3.0.0.
# Source BAM: data/benchmark/pbmc_10k_v3/pbmc_10k_v3_possorted_genome_bam.bam
# The fixture MUST keep: (a) @RG IDs containing underscores, (b) CB tags with the
# "-1" GEM-group suffix, (c) unmapped reads that still carry a CB tag,
# (d) secondary (0x100) and duplicate-flagged (0x400) reads.
set -euo pipefail
export LC_ALL=C
BAM="${1:?path to pbmc_10k_v3_possorted_genome_bam.bam}"
OUT="${2:-tests/fixtures/cellranger_pbmc_tiny.bam}"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

# 20 kb of a gene-dense, well-expressed locus; 10% subsample -> ~2 k reads.
[ -f "$BAM.bai" ] || samtools index "$BAM"   # region fetch below requires an index
samtools view -b -o "$TMP/slice.bam" "$BAM" 1:155230000-155250000
samtools view -b -s 0.10 -o "$TMP/sub.bam" "$TMP/slice.bam"

# The real BAM tags NO unmapped read with CB (500k/500k carry CR only), so the
# reference_end-is-None path must be synthesised or it is never exercised.
python3 - > "$TMP/synth.sam" <<'PY'
seq = "ACGT" * 22 + "ACG"          # 91 bp, matches this library's R2 length
for i in range(1, 21):
    print("\t".join([f"SYNTH_UNMAPPED_{i}", "4", "*", "0", "0", "*", "*", "0", "0",
                     seq, "F" * len(seq), "CB:Z:AAACCCAAGCGCCCAT-1",
                     "UB:Z:AAACCCAAGCGC", "RG:Z:pbmc_10k_v3:0:1:HFWFVDMXX:1"]))
PY

{ printf '@HD\tVN:1.4\tSO:coordinate\n'
  printf '@SQ\tSN:1\tLN:248956422\n'
  printf '@RG\tID:pbmc_10k_v3:0:1:HFWFVDMXX:1\tSM:pbmc_10k_v3\tLB:0.1\tPL:ILLUMINA\n'
  printf '@RG\tID:pbmc_10k_v3:0:1:HFWFVDMXX:2\tSM:pbmc_10k_v3\tLB:0.1\tPL:ILLUMINA\n'
  printf '@PG\tID:STAR\tPN:STAR\tVN:STAR_2.5.1b\n'
  samtools view "$TMP/sub.bam"
  cat "$TMP/synth.sam"
} > "$TMP/fixture.sam"

samtools view -b -o "$OUT" "$TMP/fixture.sam"
# .bai is NOT committed: build it in the test via pysam.index() into tmp_path.
echo "wrote $OUT ($(stat -c%s "$OUT") bytes, $(samtools view -c "$OUT") reads)"
