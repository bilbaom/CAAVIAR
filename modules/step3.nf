// modules/step3.nf
// Per-sample microhomology analysis (Python script).
//
// Input:  tuple(name, all_events_del.tsv, summary_df.tsv) from STEP2
// Output: tuple(name, <name>_merged_summary.csv)           emit: summary

process STEP3_MH {

    tag "${name}"

    publishDir "${params.outdir}/results/${params.param_dir}/${name}", mode: 'copy', overwrite: true

    input:
    tuple val(name), path(events_tsv), path(summary_tsv)
    val  amplicon               // reference amplicon sequence
    val  cut_site               // cut site position

    output:
    tuple val(name), path("${name}_merged_summary.csv"), emit: summary

    script:
    """
    step3_mh_analysis.py \\
        ${events_tsv}  \\
        '${amplicon}'  \\
        ${cut_site}    \\
        \$(pwd)

    # Rename generic output to include the sample name
    mv merged_summary.csv ${name}_merged_summary.csv
    """
}
