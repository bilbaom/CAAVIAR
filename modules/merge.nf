// modules/merge.nf
// Merge paired-end reads with FLASH

process MERGE_READS {

    tag "$name"

    publishDir "${params.outdir}/merged_reads", mode: 'copy', overwrite: true

    input:
    tuple val(name), path(r1_clean), path(r2_clean)

    output:
    tuple val(name), path("${name}.extendedFrags.fastq"), emit: merged

    script:
    """
    flash \\
        --min-overlap 10 \\
        --max-overlap 300 \\
        -t ${task.cpus - 1} \\
        -o ${name} \\
        ${r1_clean} \\
        ${r2_clean} > ${name}.flash.log 2>&1
    """
}
