// modules/map_reads.nf
// Map amplicon reads to genome with minimap2 and sort the output BAM.


process MAP_READS {

    tag "$name"

    input:
    tuple val(name), path(filtered_fastq)
    path genome_mmi

    output:
    tuple val(name), path("${name}.bam"), emit: bam

    script:
    """
    minimap2 \\
        -A ${params.minimap_A} \\
        -B ${params.minimap_B} \\
        -O ${params.minimap_O} \\
        -E ${params.minimap_E} \\
        -t ${task.cpus} \\
        -a ${genome_mmi} \\
        ${filtered_fastq} | \\
    samtools sort -o ${name}.bam
    """
}