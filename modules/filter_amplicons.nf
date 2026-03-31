// modules/filter_amplicons.nf
// Keep only reads spanning the full amplicon (strict seqkit amplicon filter)
// This removes degraded / chimeric reads that cause problems in CrispRVariants

process FILTER_AMPLICONS {

    tag "$name"

    publishDir "${params.outdir}/merged_reads", mode: 'copy', overwrite: true

    input:
    tuple val(name), path(merged_fastq)

    output:
    tuple val(name), path("${name}.extendedFrags.ont.fastq"), emit: filtered

    script:
    """
    seqkit amplicon \\
        -m 2 \\
        --strict-mode \\
        -F '${params.fw}' \\
        -R '${params.rev}' \\
        ${merged_fastq} > ${name}.extendedFrags.ont.fastq
    """
}
