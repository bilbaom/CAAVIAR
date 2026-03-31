// modules/step1_call_indels.nf
// Run CrispRVariants R script on ALL samples at once.
//
// All experiment-specific values are passed as explicit val inputs
// instead of referencing params.* inside the script block — this
// makes the module fully self-contained and reusable.
//
// Output: the results/<param_dir>/ subtree written by the R script.
//   main.nf then maps sample names onto their specific restable/instable
//   files and fans them out into a per-sample channel for STEP2.

process STEP1_CALL_INDELS {

    tag "all_samples"

    publishDir "${params.outdir}/results", mode: 'copy', overwrite: true

    input:
    path bam_files              // all BAMs staged into work dir (collected)
    path bai_files              // all BAI files staged into work dir (collected)
    val  csv                    // path to sample CSV
    val  amplicon               // reference amplicon sequence
    val  cut_site               // cut site position
    val  chrom                  // chromosome name (blat_T_name)
    val  t_start                // target region start
    val  t_end                  // target region end
    val  restable               // results output filename
    val  instable               // insertions output filename
    val  param_dir              // minimap2 param label, e.g. A5_B4_O25_E1

    output:
    // Emit the entire param_dir subtree.
    // main.nf maps sample names onto their specific files from here.
    path "results/${param_dir}",  emit: results_dir

    script:
    """
    # Recreate the bam/<param_dir>/ structure the R script expects,
    # symlinking the staged BAMs and BAIs into it.
    mkdir -p bam/${param_dir}
    for bam in *.bam; do
        [ -f "\$bam" ] || continue
        ln -sf \$(readlink -f \$bam) bam/${param_dir}/\$bam
    done
    for bai in *.bam.bai; do
        [ -f "\$bai" ] || continue
        ln -sf \$(readlink -f \$bai) bam/${param_dir}/\$bai
    done

    run_crisprvariants_universal.R \\
        \$(pwd)        \\
        ${csv}         \\
        '${amplicon}'  \\
        ${cut_site}    \\
        ${chrom}       \\
        ${t_start}     \\
        ${t_end}       \\
        ${restable}    \\
        ${instable}
    """
}
