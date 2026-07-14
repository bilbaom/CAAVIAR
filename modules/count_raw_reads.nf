// modules/count_raw_reads.nf
// Count raw reads from R1 fastq file

process COUNT_RAW_READS {
    tag { "$name" }
    publishDir { "${params.outdir}/raw_counts" }, mode: 'copy', overwrite: true

    input:
    tuple val(name), path(r1), path(r2)

    output:
    tuple val(name), path("${name}_raw_read_count.txt"), emit: count

    script:
    """
    READS=\$(zcat ${r1} | wc -l)
    echo \$(( READS / 4 )) > ${name}_raw_read_count.txt
    """
}
