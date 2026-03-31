// modules/step2.nf
// Per-sample AAV insertion analysis (R script).
//
// Input:  tuple(name, restable, instable) — the two CrispRVariants output
//         files for this sample, staged directly into the work dir.
//         No BAM file needed here.
//
// Output: tuple(name, all_events_del.tsv, summary_df.tsv)  emit: events

process STEP2_INSERTIONS {

    tag "${name}"

    publishDir "${params.outdir}/results/${params.param_dir}/${name}", mode: 'copy', overwrite: true

    input:
    tuple val(name), path(restable), path(instable)
    val  blast_db               // path to BLAST db prefix

    output:
    tuple val(name), path("all_events_del.tsv"), path("summary_df.tsv"), emit: events

    script:
    """
    export BLAST_DB_PATH="${blast_db}"

    step2_aav_insertions.R \\
        ${restable} \\
        ${instable} \\
        "${name}"
    """
}
