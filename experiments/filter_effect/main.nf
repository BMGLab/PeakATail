#!/usr/bin/env nextflow
/*
 * PeakATail — FILTER-EFFECT experiment (Laughney 6-GSM method-comparison subset)
 * ================================================================================
 * Compares how each PAS filter changes CLUSTERING, DIFFERENTIAL APA, and 3'UTR
 * SHORTEN/LENGTHEN, when the filter is applied going into the downstream matrix.
 *
 * User framing: PAS RESULTS always keep every PAS at the `ema run` annotate
 * layer -- non-matching / flagged PAS are recorded (atlas_match, internal_priming,
 * ...), never silently dropped, in the DEFAULT ("annotate") mode. Each scenario
 * below flips exactly ONE filter axis to its DROP ("filter") mode so that axis's
 * effect on the downstream matrix (clustering + switch) can be isolated,
 * one-factor-at-a-time, against the `baseline` (everything kept).
 *
 * IMPORTANT CAVEAT (see FILTER_EFFECT_README.md / the launch report): at the
 * `ema run` code level, --atlas-mode filter / --ip-filter-mode filter /
 * --annot-filter all drop the PAS BEFORE find_close() builds the gene-assignment
 * table that annotate() joins against (ema/main.py::_apply_pas_filters, called
 * before write_pas_gene_artifacts()). That means a "filter"-scenario's OWN
 * annotatedpas.bed/pasbed.bed already has fewer rows than baseline's -- the drop
 * is not confined to the clustering matrix alone. This still answers the
 * question asked (how does each filter change downstream clustering/diff/switch,
 * and how many PAS does it remove versus keeping everything) -- baseline is the
 * one scenario where the PAS list is provably complete, and each filtered
 * scenario's PAS-count delta against baseline IS the "removed for downstream"
 * number. There is no `ema` flag that keeps a PAS in the result list while
 * excluding it ONLY from clustering -- that would require new engine code,
 * out of scope here.
 *
 * SCENARIOS (all: lambda_gradient, 6-GSM subset spanning Normal/StageI/
 * IVprimary/Met -- same SUBSET_IDS as the corrected-rerun sweep's A2 tier):
 *   1. baseline           -- atlas-mode annotate, ip-filter annotate, no annot-filter
 *   2. atlas_filter        -- atlas-mode FILTER (others = baseline)
 *   3. annot_filter_3utr   -- --annot-filter --annotation-bed <3'UTR-only bed>
 *                             (others = baseline)
 *   4. ip_filter            -- --ip-filter-mode FILTER (others = baseline)
 *
 * Each scenario runs the FULL pipeline: ema run -> gex_celltyping.py ->
 * ema switch combine --group-key stage --split-key celltype -> per-celltype
 * ema switch diff (fisher, nb_multi) + switch length (classic, shannon) +
 * switch trend. Mirrors the corrected-rerun sweep's Phase-1 pattern
 * (experiments/laughney/main.nf) at 6-GSM subset scale x 4 scenarios instead
 * of 1 x 17-GSM cohort.
 *
 * DEDICATED OUTPUT ROOT (nothing under RERUN_2026-08_fixed touched):
 *   /mnt/ssd2/Laugney_Aligned/peakatail_experiments/FILTER_EFFECT_2026-08
 *
 * Run: see run_experiment.sh
 */

nextflow.enable.dsl = 2

// ── dedicated, previously-unused output root ─────────────────────────────
params.exp_root = '/mnt/ssd2/Laugney_Aligned/peakatail_experiments/FILTER_EFFECT_2026-08'
params.out_root = "${params.exp_root}/runs"
params.logs     = "${params.exp_root}/logs"

// ── tool + reference paths ────────────────────────────────────────────────
params.ema     = '/mnt/ssd2/Laugney_Aligned/.peakatail_env/bin/ema'
params.python  = '/mnt/ssd2/Laugney_Aligned/.peakatail_env/bin/python'
params.repo    = '/home/amiramiritabat/PeakATail'
// 6-GSM method-comparison subset config, REUSED read-only from the corrected
// rerun sweep's pipeline (same BAM sources / dataset ids -- never written to).
params.subset_config = '/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed/pipeline/configs/A2_base_lambda_gradient.yaml'
params.gtf     = '/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf'
params.atlas          = '/mnt/ssd2/Laugney_Aligned/refs/polyasite_3.0_GRCh38_ensembl_sorted.bed'
params.atlas_distance = 50
params.genome_fasta   = '/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.dna.primary_assembly.fa'

// GEX cell-typing inputs (same as the corrected-rerun sweep).
params.star_root = '/mnt/ssd2/Laugney_Aligned'
params.master    = '/home/amiramiritabat/laughney_realign/peakatail_laughney.yaml'
params.refs      = '/home/amiramiritabat/PeakATail/temp/laughney_refs'

// 3'UTR-only annotation BED (built once by BUILD_3UTR_BED into exp_root/refs).
params.three_utr_bed = "${params.exp_root}/refs/three_prime_utr.bed"

params.stage_order      = 'Normal,StageI,IVprimary,Met'
params.min_switch_cells = 40
params.threads_per_run  = 16
// SHORT tmp base (ssd2) -- AF_UNIX socket path-length gotcha, see corrected
// rerun sweep memory notes. Each task uses ${tmpbase}/<pid>.
params.tmpbase = '/mnt/ssd2/Laugney_Aligned/rt'

// ═════════════════════════════════════════════════════════════════════════
process PREWARM_GTF {
    tag 'gtf-cache'
    cpus 4
    output:
    val 'ready', emit: ready
    script:
    """
    export TMPDIR=${params.tmpbase}/\$\$ && mkdir -p \$TMPDIR
    ${params.python} - <<'PY'
from ema.annotate.gtf_cache import process_gtf_cached
import os
d = os.path.abspath("gtfwarm"); os.makedirs(d, exist_ok=True)
n = process_gtf_cached(gtf_path="${params.gtf}", output_dir=d,
                       endbed_path=os.path.join(d,"gene_end.bed"),
                       features_path=os.path.join(d,"raw_feature.tsv"))
print("GTF cache warmed:", len(n), "genes with UTR lengths")
PY
    """
}

// Build the 3'UTR-only annotation BED once (annot_filter_3utr scenario input).
// Extraction is a straight GTF column-3=="three_prime_utr" scan + bedtools
// merge -s (script: scripts/build_3utr_bed.py) -- independent of the GTF
// gene/isoform caches, so it can run in parallel with PREWARM_GTF.
process BUILD_3UTR_BED {
    tag '3utr-bed'
    cpus 2
    output:
    val "${params.three_utr_bed}", emit: bed
    script:
    """
    export TMPDIR=${params.tmpbase}/\$\$ && mkdir -p \$TMPDIR
    mkdir -p ${params.exp_root}/refs
    ${params.python} ${params.repo}/scripts/build_3utr_bed.py \
        --gtf ${params.gtf} --out ${params.three_utr_bed} --tmpdir \$TMPDIR
    """
}

// One `ema run` per scenario -- SAME subset config, ONE filter axis flipped.
process RUN_SCENARIO {
    tag { row.name }
    cpus { params.threads_per_run }
    input:
    tuple val(row), val(ready), val(three_utr_bed)
    output:
    tuple val(row.name), val("${params.out_root}/${row.name}"), emit: run
    script:
    def atlasflag = "--atlas ${params.atlas} --atlas-distance ${params.atlas_distance} --atlas-mode ${row.atlas_mode}"
    def ipflag    = "--ip-filter --genome-fasta ${params.genome_fasta} --ip-filter-mode ${row.ip_filter_mode}"
    def annotflag = row.annot_filter ? "--annot-filter --annotation-bed ${three_utr_bed}" : ''
    """
    export TMPDIR=${params.tmpbase}/\$\$ && mkdir -p \$TMPDIR
    ${params.ema} run -c ${params.subset_config} \
        ${atlasflag} ${ipflag} ${annotflag} \
        --output ${params.out_root}/${row.name} \
        --threads ${params.threads_per_run} --no-progress
    """
}

// GEX cell typing per scenario (same 6 GSMs each time; cheap, re-derives the
// cell-type labels each scenario's own run needs -- concordance is scenario-local).
process GEX_CELLTYPE {
    tag { name }
    cpus 16
    input:
    tuple val(name), val(run_dir)
    output:
    tuple val(name), val(run_dir), val("${run_dir}/B2_gex_celltyping"), emit: typed
    script:
    """
    export TMPDIR=${params.tmpbase}/\$\$ && mkdir -p \$TMPDIR
    ${params.python} ${params.repo}/scripts/gex_celltyping.py \
        --pas-run ${run_dir} \
        --star-root ${params.star_root} \
        --master ${params.master} \
        --refs ${params.refs} \
        --out ${run_dir}/B2_gex_celltyping
    """
}

process SWITCH_COMBINE {
    tag { name }
    cpus 8
    input:
    tuple val(name), val(run_dir), val(b2dir), val(combine_specs)
    output:
    tuple val(name), val(run_dir), path('combined/*.h5ad'), emit: celltypes
    tuple val(name), path('combined/pasbed.bed'), emit: pasbed
    script:
    """
    export TMPDIR=${params.tmpbase}/\$\$ && mkdir -p \$TMPDIR
    mkdir -p combined
    SRC=${run_dir}/annotatedpas.bed
    [ -s "\$SRC" ] || SRC=${run_dir}/pasbed.bed
    awk -F'\\t' '!s[\$4]++ {print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t0\\t"(\$7==""?\$6:\$7)}' "\$SRC" > combined/pasbed.bed
    ${params.ema} switch combine ${combine_specs} \
        --group-key stage --split-key celltype --min-cells ${params.min_switch_cells} \
        --no-progress -o combined
    """
}

process SWITCH_CELLTYPE {
    tag { "${name}/${h5.baseName}" }
    cpus 6
    input:
    tuple val(name), val(run_dir), path(h5), path(pb)
    output:
    val "${name}/${h5.baseName}", emit: done
    script:
    def sl  = h5.baseName
    def out = "${run_dir}/B3_switch"
    """
    export TMPDIR=${params.tmpbase}/\$\$ && mkdir -p \$TMPDIR
    ${params.ema} switch diff -i ${h5} --pasbed ${pb} --gtf ${params.gtf} \
        --cluster-key stage --strategy fisher --marker-top-n 0 --fdr 0.05 \
        --threads ${task.cpus} --no-progress --no-plots -o ${out}/diff/${sl}/fisher \
        || echo "WARN diff/fisher failed for ${name}/${sl}"
    ${params.ema} switch diff -i ${h5} --pasbed ${pb} --gtf ${params.gtf} \
        --cluster-key stage --strategy nb_multi --marker-top-n 200 --fdr 0.05 \
        --threads ${task.cpus} --no-progress --no-plots -o ${out}/diff/${sl}/nb_multi \
        || echo "WARN diff/nb_multi failed for ${name}/${sl}"
    for S in classic shannon; do
      ${params.ema} switch length -i ${h5} --cluster-key stage --strategy \$S \
          --pdui-pseudocount 1.0 --threads ${task.cpus} --no-progress --no-plots \
          -o ${out}/length/${sl}/\$S || echo "WARN length/\$S failed for ${name}/${sl}"
    done
    PDUI=\$(find ${out}/length/${sl}/classic -name 'pdui_classic.tsv' 2>/dev/null | head -1)
    if [ -n "\$PDUI" ]; then
      ${params.ema} switch trend --pdui "\$PDUI" \
          --stage-order ${params.stage_order} --stage-col cluster --value-col pdui \
          --no-progress -o ${out}/trend/${sl} || echo "WARN trend failed for ${name}/${sl}"
    fi
    """
}

// ═════════════════════════════════════════════════════════════════════════
workflow {
    ch_gtf_ready = PREWARM_GTF().ready.first()
    ch_3utr_bed  = BUILD_3UTR_BED().bed.first()

    // dataset_id -> group, for the 6-GSM subset (same ids/groups as the
    // corrected-rerun sweep's stage_groups.tsv, inlined here so this
    // pipeline stays self-contained -- *.tsv data tables aren't committed
    // to the repo, see .gitignore).
    id2grp = [
        'GSM3516675-Normal'         : 'Normal',
        'GSM3516663-StageIA'        : 'StageI',
        'GSM3516667-StageIA'        : 'StageI',
        'GSM3516665-StageIVprimary' : 'IVprimary',
        'GSM3516668-MetBrain'       : 'Met',
        'GSM3516664-MetBone'        : 'Met',
    ]

    scenarios = [
        [name: 'baseline',          atlas_mode: 'annotate', ip_filter_mode: 'annotate', annot_filter: false],
        [name: 'atlas_filter',      atlas_mode: 'filter',   ip_filter_mode: 'annotate', annot_filter: false],
        [name: 'annot_filter_3utr', atlas_mode: 'annotate', ip_filter_mode: 'annotate', annot_filter: true],
        [name: 'ip_filter',         atlas_mode: 'annotate', ip_filter_mode: 'filter',   annot_filter: false],
    ]

    scenario_ch = Channel.fromList(scenarios).combine(ch_gtf_ready).combine(ch_3utr_bed)
    runs  = RUN_SCENARIO(scenario_ch).run
    typed = GEX_CELLTYPE(runs)

    combine_ch = typed.map { name, run_dir, b2dir ->
        def specs = id2grp.collect { gsm, group -> "-i ${group}=${b2dir}/${gsm}_pas_labeled.h5ad" }.join(' ')
        tuple(name, run_dir, b2dir, specs)
    }
    combined = SWITCH_COMBINE(combine_ch)
    ct = combined.celltypes
            .flatMap { name, run_dir, files ->
                (files instanceof List ? files : [files])
                    .findAll { it.name.endsWith('.h5ad') }
                    .collect { tuple(name, run_dir, it) }
            }
            .combine(combined.pasbed, by: 0)   // join on scenario name -> that scenario's pasbed
            .map { name, run_dir, h5, pb -> tuple(name, run_dir, h5, pb) }
    SWITCH_CELLTYPE(ct)
}
