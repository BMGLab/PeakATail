#!/usr/bin/env nextflow
/*
 * PeakATail — Laughney 2020 LUAD cohort sweep (CORRECTED RERUN)
 * ============================================================
 * End-to-end, package-command-only orchestration of the Laughney re-analysis on
 * the 112-core BioLab node. Every heavy step is a real `ema` subcommand — no
 * pipeline logic lives here, only sequencing + fan-out + isolation.
 *
 * WHY THIS EXISTS
 * ---------------
 * All prior Laughney results were produced by buggy May code (≈50% cell
 * doubling in concat_matrices, atlas HARD-snap dropping 43–68% of PAS, E4 ledger
 * mis-resolution, categorical-key crash in switch diff). Those are fixed on
 * `develop` (PR #58). This sweep regenerates the *trustworthy* results from
 * scratch into a DEDICATED, previously-unused output root — nothing under the
 * old runs/ or runs_sweep/ is read or written.
 *
 * DEDICATED OUTPUT ROOT (params.sweep_root):
 *   /mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed
 * Everything (peak-calls, GEX typing, switches, strategy runs, reannotate
 * branches, benchmarks, nextflow logs) nests inside it. See README + run_sweep.sh.
 *
 * TWO PHASES (get the headline result EARLY, then fan out)
 * -------------------------------------------------------
 *  PHASE 1  (critical path — the professor-ready headline):
 *    lambda_gradient, FULL 17-GSM cohort
 *      → ema run                       (peak-call + cluster, atlas ANNOTATE, no drop)
 *      → gex_celltyping.py             (GEX cell types per GSM; parked A1 script)
 *      → ema switch combine            (per-celltype, stage-labelled h5ads)
 *      → ema switch diff/length/trend  (APA shift + 3'UTR shortening ACROSS STAGES,
 *                                       within each cell type — the biological Q)
 *      → ema switch match              (cross-dataset cell-type sanity check)
 *
 *  PHASE 2  (runs in parallel once the cohort peak-call exists):
 *    other strategies + one-factor-at-a-time filter/cluster sweep
 *      → ema run  (lambda_poisson, sierra_iterative)  → bench vs PolyASite v3.0
 *      → ema reannotate  (trim/cluster branches off the cohort's SAVED peak-calls
 *                         — never re-peak-calls)
 *
 * Select with --phase {phase1|phase2|all}. Default 'all' gates phase-2 reannotate
 * on the cohort peak-call so the two overlap without racing the same files.
 *
 * ANNOTATE-NOT-DROP (user requirement):
 *   atlas + internal-priming are OVERLAYS. `--atlas-mode annotate` /
 *   `--ip-filter-mode annotate` KEEP every PAS and stamp status columns
 *   (atlas_match / atlas_distance_bp / internal_priming). Sweeps therefore never
 *   discard non-PolyASite ("alternative polyA") sites. Match/no-match counts are
 *   emitted at every run (atlas_stats.json) and every benchmark.
 *
 * ISOLATION (unchanged, still holds):
 *   1. Each `ema run`/branch writes its OWN absolute output_dir — no shared writes.
 *   2. ema's process-global singletons live inside each task's own `ema`
 *      subprocess — no cross-task state leak.
 *   3. The one shared-WRITE resource (global GTF cache) is populated ONCE by
 *      PREWARM_GTF before any fan-out; every later task only READS it.
 *   4. Every task exports a private TMPDIR (pybedtools/samtools scratch).
 *   5. Barcodes namespaced by hyphenated dataset id — verified collision-free.
 *
 * Run:  see run_sweep.sh   (nextflow run main.nf -profile biolab ...)
 */

nextflow.enable.dsl = 2

// ── dedicated, previously-unused output root ─────────────────────────────
params.sweep_root = '/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed'
params.out_root   = "${params.sweep_root}/runs"
params.logs       = "${params.sweep_root}/logs"

// ── tool + reference paths (override on CLI or in nextflow.config) ────────
params.ema     = '/mnt/ssd2/Laugney_Aligned/.peakatail_env/bin/ema'
params.python  = '/mnt/ssd2/Laugney_Aligned/.peakatail_env/bin/python'
params.repo    = '/home/amiramiritabat/PeakATail'
params.configs = "${projectDir}/configs"
params.gtf     = '/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf'

// PolyASite v3.0 (2025, scRNA-derived, GRCh38, Ensembl-reconciled + sorted).
// Used ONLY as an annotate overlay + benchmark reference — never a drop filter.
params.atlas          = '/mnt/ssd2/Laugney_Aligned/refs/polyasite_3.0_GRCh38_ensembl_sorted.bed'
params.atlas_distance = 50
// Internal-priming (fasta) annotate overlay. Phase 1 keeps this ON in ANNOTATE
// mode — every PAS near an A-rich stretch is KEPT and flagged (internal_priming),
// used for clustering + switch like any other PAS. Override to '' to disable.
params.genome_fasta   = '/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.dna.primary_assembly.fa'

// GEX cell-typing inputs (gex_celltyping.py — the parked A1 GEX-I/O script).
params.star_root = '/mnt/ssd2/Laugney_Aligned'   // holds <SRR>_STAR/<SRR>_Solo.out/Gene/filtered
params.master    = '/home/amiramiritabat/laughney_realign/peakatail_laughney.yaml'
params.refs      = '/home/amiramiritabat/PeakATail/temp/laughney_refs'  // marker-signature CSVs (read-only)

// experiment tables
params.groups        = "${projectDir}/stage_groups.tsv"        // GSM → group,stage
params.grid          = "${projectDir}/trim_cluster_grid.tsv"   // reannotate branches
params.peakcall_grid = "${projectDir}/peakcall_grid.tsv"       // ema-run method sweep

// biological ordering for the 3'UTR-length stage trend (well-powered groups).
// SHORT tmp base (ssd2). MUST be short: ema's multiprocessing Manager opens an
// AF_UNIX socket at $TMPDIR/pymp-*/listener-* and the deep work-dir path overflows
// the 107-char AF_UNIX limit. Each task uses ${tmpbase}/<pid> (unique + short).
params.tmpbase          = '/mnt/ssd2/Laugney_Aligned/rt'
params.stage_order      = 'Normal,StageI,IVprimary,Met'
params.min_switch_cells = 40      // drop a (celltype,stage) cell with fewer cells

// ── PHASE 3 (FILTER-EFFECT experiment) params — dedicated output root,
// nothing under params.sweep_root/params.out_root is touched by phase3. ──
params.filter_effect_root     = '/mnt/ssd2/Laugney_Aligned/peakatail_experiments/FILTER_EFFECT_2026-08'
params.filter_effect_out      = "${params.filter_effect_root}/runs"
// The subset base run phase3 branches off (produced by phase2's GRID_RUN,
// row "lg_annotate" in peakcall_grid.tsv: 6-GSM subset, lambda_gradient,
// atlas annotate + ip annotate — "keep everything"). READ-ONLY.
params.filter_effect_base_run = "${params.out_root}/grid/lg_annotate"
params.three_utr_bed          = "${params.filter_effect_root}/refs/three_prime_utr.bed"

// phase selection (phase-2 method axes live in peakcall_grid.tsv + trim_cluster_grid.tsv)
params.phase = 'all'   // phase1 | phase2 | phase3 | all  (phase3 is opt-in ONLY — see below)

def do_p1 = (params.phase == 'all' || params.phase == 'phase1')
def do_p2 = (params.phase == 'all' || params.phase == 'phase2')
// PHASE 3 — FILTER-EFFECT experiment. Deliberately EXCLUDED from 'all' (opt-in
// only, `--phase phase3`) so a routine `bash run_sweep.sh -resume` NEVER
// triggers it by surprise. Continues from phase2's already-computed 6-GSM
// subset peak-calls (runs/grid/lg_annotate — lambda_gradient, atlas
// annotate, ip annotate) via `ema reannotate` — NO re-peak-calling.
//
// Every scenario LABELS all three axes (atlas_match, internal_priming,
// in_3utr) on every PAS -- annotatedpas.bed carries the FULL labeled PAS
// set, identically, across all 4 scenarios (genuinely "keeps everything",
// not just lg_annotate's own copy). Independently, each scenario flips
// EXACTLY ONE `--exclude-*` flag (or none, for baseline), which narrows
// ONLY the clustering matrix that scenario's own GEX-celltyping/diff/
// length/trend re-run on -- see FILTER_EFFECT_BRANCH + ema/reannotate.py's
// "PAS labels + clustering mask" for the full mechanism. lg_annotate's own
// PAS results (posbed.bed/negbed.bed/annotatedpas.bed) are never touched.
// Output goes to a DEDICATED root (params.filter_effect_root) — nothing
// under params.out_root/params.sweep_root is written by phase3.
def do_p3 = (params.phase == 'phase3')

// ═════════════════════════════════════════════════════════════════════════
// PHASE 0 — prewarm the global GTF cache ONCE (serial) so peak-calls only read
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

// ═════════════════════════════════════════════════════════════════════════
// PHASE 1 — lambda_gradient, full cohort, headline result
// ═════════════════════════════════════════════════════════════════════════

// 1.1 the one expensive peak-call: 17 GSMs, lambda_gradient, atlas ANNOTATE.
process RUN_COHORT {
    tag 'cohort-lambda_gradient'
    cpus { params.threads_per_run }
    input:
    val ready
    output:
    val "${params.out_root}/B1_cohort_full", emit: cohort
    script:
    // atlas + fasta BOTH annotate-mode: keep every non-matching PAS, use all for
    // clustering + switch. --ip-filter ENABLES the fasta step; annotate keeps all.
    def ip = params.genome_fasta ? "--ip-filter --genome-fasta ${params.genome_fasta} --ip-filter-mode annotate" : ''
    """
    export TMPDIR=${params.tmpbase}/\$\$ && mkdir -p \$TMPDIR
    ${params.ema} run -c ${params.configs}/B1_cohort_full.yaml \
        --atlas ${params.atlas} --atlas-distance ${params.atlas_distance} --atlas-mode annotate ${ip} \
        --output ${params.out_root}/B1_cohort_full \
        --threads ${params.threads_per_run} --no-progress
    """
}

// 1.2 GEX cell typing per GSM (pools STARsolo GEX, marker-scores cell types,
//     writes <gsm>_pas_labeled.h5ad with obs['celltype'] + concordance stats).
process GEX_CELLTYPE {
    tag 'gex-celltyping'
    cpus 16
    input:
    val cohort
    output:
    tuple val(cohort), val("${cohort}/B2_gex_celltyping"), emit: typed
    script:
    """
    export TMPDIR=${params.tmpbase}/\$\$ && mkdir -p \$TMPDIR
    ${params.python} ${params.repo}/scripts/gex_celltyping.py \
        --pas-run ${cohort} \
        --star-root ${params.star_root} \
        --master ${params.master} \
        --refs ${params.refs} \
        --out ${cohort}/B2_gex_celltyping
    """
}

// 1.3 combine per-GSM labelled h5ads into per-CELLTYPE, stage-labelled h5ads.
//     A complete pasbed.bed (from annotatedpas.bed — the incomplete root pasbed
//     has only ~5.8k of ~30k PAS) is placed as a sibling so switch length/diff
//     resolve full coords by walking up from each combined h5ad.
process SWITCH_COMBINE {
    tag 'switch-combine'
    cpus 8
    input:
    tuple val(cohort), val(b2dir), val(combine_specs)
    output:
    tuple val(cohort), path('combined/*.h5ad'), emit: celltypes
    path 'combined/pasbed.bed', emit: pasbed
    script:
    """
    export TMPDIR=${params.tmpbase}/\$\$ && mkdir -p \$TMPDIR
    mkdir -p combined
    # complete pasbed sibling (chrom,start,end,pas_id,score,strand); dedup pas_id.
    SRC=${cohort}/annotatedpas.bed
    [ -s "\$SRC" ] || SRC=${cohort}/pasbed.bed
    awk -F'\\t' '!s[\$4]++ {print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t0\\t"(\$7==""?\$6:\$7)}' "\$SRC" > combined/pasbed.bed
    ${params.ema} switch combine ${combine_specs} \
        --group-key stage --split-key celltype --min-cells ${params.min_switch_cells} \
        --no-progress -o combined
    """
}

// 1.4 per cell type: differential APA (fisher + nb_multi omnibus across stages),
//     3'UTR length (classic PDUI + proportion + shannon), ordered-stage trend.
process SWITCH_CELLTYPE {
    tag { h5.baseName }
    cpus 6
    input:
    tuple val(cohort), path(h5), path(pb)
    output:
    val "${h5.baseName}", emit: done
    script:
    def sl  = h5.baseName
    def out = "${cohort}/B3_switch"
    // pb (the complete pasbed) is STAGED next to h5 in this work dir, so diff reads
    // it via --pasbed and length resolves it by walking up from the h5ad sibling.
    // Each analysis is fault-tolerant: a failing strategy on one cell type must not
    // drop the others (WARN is visible in the log, not silently swallowed).
    """
    export TMPDIR=${params.tmpbase}/\$\$ && mkdir -p \$TMPDIR
    # diff — fisher = exhaustive within-gene screen; nb_multi = omnibus LRT over all stages
    ${params.ema} switch diff -i ${h5} --pasbed ${pb} --gtf ${params.gtf} \
        --cluster-key stage --strategy fisher --marker-top-n 0 --fdr 0.05 \
        --threads ${task.cpus} --no-progress --no-plots -o ${out}/diff/${sl}/fisher \
        || echo "WARN diff/fisher failed for ${sl}"
    ${params.ema} switch diff -i ${h5} --pasbed ${pb} --gtf ${params.gtf} \
        --cluster-key stage --strategy nb_multi --marker-top-n 200 --fdr 0.05 \
        --threads ${task.cpus} --no-progress --no-plots -o ${out}/diff/${sl}/nb_multi \
        || echo "WARN diff/nb_multi failed for ${sl}"
    # length — PDUI (classic) + full usage vector (proportion) + entropy (shannon)
    for S in classic proportion shannon; do
      ${params.ema} switch length -i ${h5} --cluster-key stage --strategy \$S \
          --pdui-pseudocount 1.0 --threads ${task.cpus} --no-progress --no-plots \
          -o ${out}/length/${sl}/\$S || echo "WARN length/\$S failed for ${sl}"
    done
    # trend — ordered-stage 3'UTR shortening/lengthening slope + Spearman
    PDUI=\$(find ${out}/length/${sl}/classic -name 'pdui_classic.tsv' 2>/dev/null | head -1)
    if [ -n "\$PDUI" ]; then
      ${params.ema} switch trend --pdui "\$PDUI" \
          --stage-order ${params.stage_order} --stage-col cluster --value-col pdui \
          --no-progress -o ${out}/trend/${sl} || echo "WARN trend failed for ${sl}"
    fi
    """
}

// 1.5 cross-dataset cell-type matching (native package check vs the GEX typing).
process SWITCH_MATCH {
    tag 'cohort-match'
    cpus 12
    input:
    val cohort
    output:
    val 'matched', emit: done
    script:
    """
    export TMPDIR=${params.tmpbase}/\$\$ && mkdir -p \$TMPDIR
    H5ADS=\$(ls ${cohort}/07_clustering/*/clusters.h5ad)
    ${params.ema} switch match \$(for h in \$H5ADS; do echo -n "-i \$h "; done) \
        --strategy marker_overlap --no-progress --no-plots \
        -o ${cohort}/B3_switch/match
    """
}

// ═════════════════════════════════════════════════════════════════════════
// PHASE 2 — other strategies + reannotate OFAT sweep (parallel)
// ═════════════════════════════════════════════════════════════════════════

// 2.1 GRID_RUN — the peak-call method sweep. Each row of peakcall_grid.tsv is one
//     `ema run` on the 6-GSM method-comparison subset, varying the method/param at
//     a step: peak_strategy × atlas_mode × ip_mode. All off ONE subset base config
//     (A2_base) via CLI overrides + a per-row --output, so nothing collides. The
//     fasta axis compares ANNOTATE (keep+flag) vs FILTER (drop internal-priming,
//     don't cluster them) vs OFF. Atlas stays ANNOTATE everywhere — atlas OVERLAP
//     is read from each run's atlas_stats.json + the vs-v3 benchmark, so no
//     atlas-filter clustering variant is run.
process GRID_RUN {
    tag { row.run_name }
    cpus { params.threads_per_run }
    input:
    tuple val(row), val(ready)
    output:
    val "${params.out_root}/grid/${row.run_name}", emit: run
    script:
    def atlasflag = (row.atlas_mode == 'off') ? '' :
        "--atlas ${params.atlas} --atlas-distance ${params.atlas_distance} --atlas-mode ${row.atlas_mode}"
    def ipflag = (row.ip_mode == 'off' || !params.genome_fasta) ? '' :
        "--ip-filter --genome-fasta ${params.genome_fasta} --ip-filter-mode ${row.ip_mode}"
    """
    export TMPDIR=${params.tmpbase}/\$\$ && mkdir -p \$TMPDIR
    ${params.ema} run -c ${params.configs}/A2_base_lambda_gradient.yaml \
        --peak-strategy ${row.peak_strategy} ${atlasflag} ${ipflag} \
        --output ${params.out_root}/grid/${row.run_name} \
        --threads ${params.threads_per_run} --no-progress
    """
}

// 2.2 benchmark each strategy's called PAS vs PolyASite v3.0 (annotate stats,
//     precision/recall at distance cutoffs — NOT a filter).
process BENCH_STRATEGY {
    tag { run.tokenize('/').last() }
    cpus 4
    input:
    val run
    output:
    val 'benched', emit: done
    script:
    def name = run.tokenize('/').last()
    """
    export TMPDIR=${params.tmpbase}/\$\$ && mkdir -p \$TMPDIR
    ${params.python} ${params.repo}/scripts/bench_pas_vs_atlas.py \
        --pasbed ${run}/pasbed.bed --atlas ${params.atlas} \
        --out ${run}/benchmark_vs_polyasite_v3.json --label ${name}
    """
}

// 2.3 trim/cluster branches off the cohort's SAVED peak-calls (ema reannotate —
//     never re-peak-calls). One branch per trim_cluster_grid.tsv row.
process REANNOTATE_BRANCH {
    tag { row.branch_name }
    cpus 8
    input:
    tuple val(row), val(cohort)
    output:
    val row.branch_name, emit: done
    script:
    def ext = row.include_extended?.toString()?.toLowerCase() == 'true' ? '--include-extended' : ''
    def nn  = row.n_neighbors ? "--n-neighbors ${row.n_neighbors}" : ''
    """
    export TMPDIR=${params.tmpbase}/\$\$ && mkdir -p \$TMPDIR
    ${params.ema} reannotate \
        --base-run ${cohort} \
        --out ${params.out_root}/reannotate/${row.branch_name} \
        --gtf ${params.gtf} \
        --max-gene-distance ${row.max_gene_distance} \
        --utr-multiplier ${row.utr_multiplier} ${ext} \
        --cluster-method ${row.cluster_method} \
        --resolution ${row.resolution} ${nn} \
        --threads ${task.cpus}
    """
}

// ═════════════════════════════════════════════════════════════════════════
// PHASE 3 — FILTER-EFFECT experiment (opt-in, `--phase phase3` only)
// ═════════════════════════════════════════════════════════════════════════

// 3.1 3'UTR-only annotation BED for the annot_filter_3utr scenario — the
//     default `ema reannotate --annot-filter` fallback (no --annotation-bed)
//     is a whole-gene-body BED (ema/annotate/gtftobed.py, source_type="gene"),
//     not 3'UTR-specific, so this is built explicitly from the GTF's
//     three_prime_utr feature rows (scripts/build_3utr_bed.py).
process PREP_3UTR_BED {
    tag '3utr-bed'
    cpus 2
    output:
    val "${params.three_utr_bed}", emit: bed
    script:
    """
    export TMPDIR=${params.tmpbase}/\$\$ && mkdir -p \$TMPDIR
    mkdir -p ${params.filter_effect_root}/refs
    ${params.python} ${params.repo}/scripts/build_3utr_bed.py \
        --gtf ${params.gtf} --out ${params.three_utr_bed} --tmpdir \$TMPDIR
    """
}

// 3.2 One `ema reannotate` branch per scenario, off the SAME saved
//     lg_annotate peak-calls (params.filter_effect_base_run) — NEVER
//     re-peak-calls. EVERY branch labels ALL THREE axes (atlas_match,
//     internal_priming, in_3utr) via --atlas/--genome-fasta/--annotation-bed
//     -- annotatedpas.bed carries the FULL labeled PAS set identically
//     across all 4 scenarios (the "keeps everything" contract). `baseline`
//     sets no --exclude-* flag (reproduces lg_annotate's own downstream
//     byte-for-byte); the other 3 each flip exactly ONE exclude flag, which
//     narrows ONLY the clustering matrix that scenario's GEX/diff/length/
//     trend re-run on. Every branch's `--out` is under
//     params.filter_effect_out — lg_annotate's own dir
//     (params.filter_effect_base_run) is read-only, never written.
process FILTER_EFFECT_BRANCH {
    tag { row.name }
    cpus 8
    input:
    tuple val(row), val(base_run), val(three_utr_bed)
    output:
    val "${params.filter_effect_out}/${row.name}", emit: run
    script:
    def labelFlags = "--atlas ${params.atlas} --atlas-distance ${params.atlas_distance} " +
                     "--genome-fasta ${params.genome_fasta} --annotation-bed ${three_utr_bed}"
    def excludeFlag = row.exclude ? "--exclude-${row.exclude}" : ''
    """
    export TMPDIR=${params.tmpbase}/\$\$ && mkdir -p \$TMPDIR
    ${params.ema} reannotate \
        --base-run ${base_run} \
        --out ${params.filter_effect_out}/${row.name} \
        --gtf ${params.gtf} \
        ${labelFlags} ${excludeFlag} \
        --threads ${task.cpus}
    """
}

// ═════════════════════════════════════════════════════════════════════════
workflow {
    ch_ready = PREWARM_GTF().ready.first()

    // group -> [dataset_ids]  and  dataset_id -> group   (from stage_groups.tsv)
    grp2ids = [:].withDefault { [] }
    id2grp  = [:]
    file(params.groups).eachLine { ln ->
        if( ln && !ln.startsWith('#') && !ln.startsWith('dataset_id') ) {
            def f = ln.split('\t'); grp2ids[f[1]] << f[0]; id2grp[f[0]] = f[1]
        }
    }

    cohort_ch = Channel.empty()

    if( do_p1 ) {
        cohort_ch = RUN_COHORT(ch_ready).cohort
        typed     = GEX_CELLTYPE(cohort_ch)

        // build combine specs: one  -i <group>=<b2>/<gsm>_pas_labeled.h5ad  per GSM
        combine_ch = typed.map { cohort, b2dir ->
            def specs = id2grp.collect { gsm, group ->
                "-i ${group}=${b2dir}/${gsm}_pas_labeled.h5ad"
            }.join(' ')
            tuple(cohort, b2dir, specs)
        }
        combined = SWITCH_COMBINE(combine_ch)
        ct = combined.celltypes
                .flatMap { cohort, files ->
                    (files instanceof List ? files : [files])
                        .findAll { it.name.endsWith('.h5ad') }
                        .collect { tuple(cohort, it) }
                }
                .combine(combined.pasbed)   // append the one complete pasbed to each celltype
        SWITCH_CELLTYPE(ct)
        SWITCH_MATCH(cohort_ch)
    }

    if( do_p2 ) {
        // peak-call method sweep (strategy × atlas_mode × ip_mode) from the grid.
        grid_run_ch = Channel.fromPath(params.peakcall_grid)
                             .splitCsv(header:true, sep:'\t')
                             .filter { !it.run_name.startsWith('#') }
                             .combine(ch_ready)
        grid_runs = GRID_RUN(grid_run_ch).run
        BENCH_STRATEGY(grid_runs)

        // reannotate branches need the cohort peak-calls. In 'all' we gate on the
        // live cohort signal; in 'phase2' the cohort must already exist on disk.
        base = do_p1 ? cohort_ch.first()
                     : Channel.value("${params.out_root}/B1_cohort_full")
        grid = Channel.fromPath(params.grid)
                      .splitCsv(header:true, sep:'\t')
                      .filter { !it.branch_name.startsWith('#') }
        REANNOTATE_BRANCH( grid.combine(base) )
    }

    if( do_p3 ) {
        // lg_annotate must already exist on disk (phase2, run separately —
        // phase3 is deliberately its own `--phase phase3` invocation so it
        // never races/duplicates a live phase2 grid run). Same 6-GSM subset
        // as gen_configs.py's SUBSET_IDS, restricted from the full id2grp.
        subset_ids = [
            'GSM3516675-Normal', 'GSM3516663-StageIA', 'GSM3516667-StageIA',
            'GSM3516665-StageIVprimary', 'GSM3516668-MetBrain', 'GSM3516664-MetBone',
        ]
        subset_id2grp = id2grp.subMap(subset_ids)

        // `exclude` maps to `--exclude-<value>` on ema reannotate (empty = no
        // exclusion, i.e. baseline). Every scenario labels all 3 axes
        // regardless (see FILTER_EFFECT_BRANCH's labelFlags) -- only which
        // ONE label (if any) also excludes from clustering differs.
        scenarios = [
            [name: 'baseline',          exclude: ''],
            [name: 'atlas_filter',      exclude: 'atlas-nonmatch'],
            [name: 'annot_filter_3utr', exclude: 'not-in-3utr'],
            [name: 'ip_filter',         exclude: 'internal-priming'],
        ]

        // FILTER_EFFECT_BRANCH/GEX_CELLTYPE/SWITCH_COMBINE reuse the exact
        // same `val cohort` / tuple shapes phase1 uses for B1_cohort_full —
        // here `cohort` is one of the 4 scenario run dirs (its basename IS
        // the scenario name, by construction of --out above), so all three
        // processes are reused UNCHANGED; only the pasbed join needs `by: 0`
        // (phase1 only ever has one cohort in flight, so a bare `.combine()`
        // broadcast was fine there — phase3 has 4 concurrent cohorts, so the
        // join MUST be keyed or celltypes from one scenario would pair with
        // another scenario's pasbed).
        ch_3utr_bed = PREP_3UTR_BED().bed.first()
        base_run_ch = Channel.value(params.filter_effect_base_run)
        fe_ch = Channel.fromList(scenarios).combine(base_run_ch).combine(ch_3utr_bed)
        fe_runs  = FILTER_EFFECT_BRANCH(fe_ch).run
        fe_typed = GEX_CELLTYPE(fe_runs)

        fe_combine_ch = fe_typed.map { cohort, b2dir ->
            def specs = subset_id2grp.collect { gsm, group ->
                "-i ${group}=${b2dir}/${gsm}_pas_labeled.h5ad"
            }.join(' ')
            tuple(cohort, b2dir, specs)
        }
        fe_combined = SWITCH_COMBINE(fe_combine_ch)
        fe_ct = fe_combined.celltypes
                .flatMap { cohort, files ->
                    (files instanceof List ? files : [files])
                        .findAll { it.name.endsWith('.h5ad') }
                        .collect { tuple(cohort, it) }
                }
                .combine(fe_combined.pasbed, by: 0)
        SWITCH_CELLTYPE(fe_ct)
    }
}
