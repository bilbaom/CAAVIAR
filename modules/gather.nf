// modules/gather.nf
// Collect all per-sample merged_summary.csv files and _aavcount.txt files
// into a single all_results_merged_summary.csv.
//
// Both input collections are passed explicitly so Nextflow waits for every
// upstream job before running gather_results.py.  The script itself locates
// the files on disk via the outdir path — the inputs here are only used to
// establish the dependency edges.

process GATHER_RESULTS {

    tag "gather"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path summaries    // all merged_summary.csv files (collected)
    path aav_counts   // all *_aavcount.txt files (collected)

    output:
    path "all_results_merged_summary.csv", emit: gathered_csv

    script:
    """
    gather_results.py ${params.outdir}
    """
}
