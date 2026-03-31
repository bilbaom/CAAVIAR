// modules/index_genome.nf
// Create a minimap2 index (.mmi) from a reference genome FASTA

process INDEX_GENOME {

    tag "${genome_fasta.baseName}"

    // Optionally publish the index so you can reuse it in future runs without rebuilding
    publishDir "${params.outdir}/reference", mode: 'copy', overwrite: true

    input:
    path genome_fasta

    output:
    path "${genome_fasta.baseName}.mmi", emit: mmi

    script:
    """
    minimap2 \\
        -t ${task.cpus} \\
        -d ${genome_fasta.baseName}.mmi \\
        ${genome_fasta}
    """
}