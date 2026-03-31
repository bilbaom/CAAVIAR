// modules/trim.nf
// Adapter and quality trimming with bbduk

process TRIM_READS {

    tag "$name"   // shows sample name in the log  e.g. [TRIM_READS (SRR123)]

    publishDir "${params.outdir}/bbtools_cleaned", mode: 'copy', overwrite: true

    input:
    tuple val(name), path(r1), path(r2)

    output:
    tuple val(name), path("${name}_clean_R1.fq.gz"), path("${name}_clean_R2.fq.gz"), emit: trimmed

    script:
    """
    bbduk.sh \\
        in1=${r1} \\
        in2=${r2} \\
        out1=${name}_clean_R1.fq.gz \\
        out2=${name}_clean_R2.fq.gz \\
        ref=adapters \\
        ktrim=r k=23 mink=11 hdist=1 \\
        qtrim=rl trimq=20 minlen=50 \\
        tpe tbo \\
        threads=${task.cpus}
    """
}
