// modules/filter_mappings.nf
// Filter a sorted BAM in two steps:
//   1. Remove supplementary (0x800) and secondary (0x100) alignments
//   2. Keep only on-target reads using a two-window samtools view
//      (reads must overlap both the start and end of the target region)

process FILTER_MAPPINGS {

    tag "$name"

    publishDir "${params.outdir}/bam/${params.param_dir}", mode: 'copy', overwrite: true

    input:
    tuple val(name), path(bam)

    output:
    tuple val(name), path("${name}_mapped_clean.bam"),     emit: bam_clean
	tuple val(name), path("${name}_mapped_clean.bam.bai"), emit: bam_index

    script:
    def window_start = params.blat_T_start.toLong()
    def window_end   = params.blat_T_end.toLong()
    """
    # 1. Remove supplementary (0x800) and secondary (0x100) alignments
    samtools view -b -F 2304 ${bam} > ${name}_clean.bam

    # 2. Filter to on-target region (two-window approach):
    #    first pass keeps reads overlapping the start window,
    #    second pass keeps only those also overlapping the end window
    samtools index ${name}_clean.bam

    samtools view -b -h ${name}_clean.bam \\
        ${params.blat_T_name}:${window_start}-\$((${window_start} + 20)) \\
        > ${name}_tmp.bam
    samtools index ${name}_tmp.bam

    samtools view -b -h ${name}_tmp.bam \\
        ${params.blat_T_name}:\$((${window_end} - 20))-${window_end} \\
        > ${name}_mapped_clean.bam
    samtools index ${name}_mapped_clean.bam
    """
}
