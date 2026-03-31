// modules/blast_aav.nf
// Convert merged FASTQ → FASTA, BLAST against AAV sequences,
// count unique hits → write to *_aavcount.txt

process BLAST_AAV {

    tag "$name"

    publishDir "${params.outdir}/merged_reads", mode: 'copy', overwrite: true

    input:
    tuple val(name), path(filtered_fastq)

    output:
    tuple val(name), path("${name}_aavcount.txt"), emit: aav_count

    script:
    """
    seqkit fq2fa ${filtered_fastq} | \\
    blastn \\
        -db ${params.blast_db} \\
        -task blastn \\
        -word_size 12 \\
        -evalue 0.01 \\
        -perc_identity 80 \\
        -max_target_seqs 10000000 \\
        -outfmt 6 \\
        -num_threads ${task.cpus} | \\
    cut -f 1 | sort -u | wc -l > ${name}_aavcount.txt
    """
}
