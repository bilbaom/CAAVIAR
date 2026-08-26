// modules/report.nf
// Generate HTML report with matplotlib plots

process GENERATE_REPORT {
    tag { "report" }
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path summary_csv
    path alignment_pngs
    path params_json
    path samples_csv

    output:
    path "pipeline_report.html", emit: report_html

    script:
    """
    generate_html_report.py ${summary_csv} pipeline_report.html ${params_json} ${samples_csv}
    """
}

