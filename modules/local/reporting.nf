process BUILD_QC_DECISION {
    tag "${sample}"
    publishDir "${params.outdir}/tables", mode: 'copy', pattern: '*.sample_summary.*'
    publishDir "${params.outdir}/tables", mode: 'copy', pattern: '*.tumor_purity_estimate.tsv'
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*.qc_flags.tsv'
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*.classification.tsv'
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*.executive_summary.md'
    publishDir "${params.outdir}/tables", mode: 'copy', pattern: '*.panel_variant_status.tsv'

    input:
    tuple val(sample), path(bam_metrics), path(coverage_summary), path(gene_coverage), path(interval_coverage), path(vcf_summary), path(variants), path(tumor_purity), path(panel_genes)
    path rules_json

    output:
    tuple val(sample), path("${sample}.sample_summary.tsv"), emit: sample_summary_tsv
    tuple val(sample), path("${sample}.sample_summary.csv"), emit: sample_summary_csv
    tuple val(sample), path("${sample}.qc_flags.tsv"), emit: qc_flags
    tuple val(sample), path("${sample}.classification.tsv"), emit: classification
    tuple val(sample), path("${sample}.executive_summary.md"), emit: executive_summary
    tuple val(sample), path("${sample}.panel_variant_status.tsv"), emit: panel_variant_status
    tuple val(sample), path("${sample}.tumor_purity_estimate.tsv"), emit: tumor_purity

    script:
    """
    ${params.python_executable} ${projectDir}/bin/build_sample_qc.py \
        --sample ${sample} \
        --bam-metrics ${bam_metrics} \
        --coverage-summary ${coverage_summary} \
        --gene-coverage ${gene_coverage} \
        --interval-coverage ${interval_coverage} \
        --vcf-summary ${vcf_summary} \
        --variants ${variants} \
        --tumor-purity ${tumor_purity} \
        --panel-genes ${panel_genes} \
        --rules ${rules_json} \
        --summary-tsv ${sample}.sample_summary.tsv \
        --summary-csv ${sample}.sample_summary.csv \
        --flags-out ${sample}.qc_flags.tsv \
        --classification-out ${sample}.classification.tsv \
        --executive-out ${sample}.executive_summary.md \
        --panel-status-out ${sample}.panel_variant_status.tsv \
        --tumor-purity-out ${sample}.tumor_purity_estimate.tsv
    """
}

process RENDER_REPORT {
    tag "${sample}"
    publishDir "${params.outdir}/reports", mode: 'copy', pattern: '*.qc_report.html'
    publishDir "${params.outdir}/reports", mode: 'copy', pattern: '*.qc_report.md'
    publishDir "${params.outdir}/reports", mode: 'copy', pattern: '*.plot_manifest.tsv'

    input:
    tuple val(sample), path(sample_summary), path(classification), path(flags), path(gene_coverage), path(interval_coverage), path(variants), path(mapq_distribution), path(tumor_purity), path(panel_variant_status), path(panel_metadata)

    output:
    tuple val(sample), path("${sample}.qc_report.html"), emit: report_html
    tuple val(sample), path("${sample}.qc_report.md"), emit: report_md
    tuple val(sample), path("${sample}.plot_manifest.tsv"), emit: plot_manifest

    script:
    """
    ${params.python_executable} ${projectDir}/bin/render_report.py \
        --sample ${sample} \
        --sample-summary ${sample_summary} \
        --classification ${classification} \
        --qc-flags ${flags} \
        --gene-coverage ${gene_coverage} \
        --interval-coverage ${interval_coverage} \
        --variants ${variants} \
        --mapq-distribution ${mapq_distribution} \
        --tumor-purity ${tumor_purity} \
        --panel-variant-status ${panel_variant_status} \
        --panel-metadata ${panel_metadata} \
        --html-out ${sample}.qc_report.html \
        --markdown-out ${sample}.qc_report.md \
        --plot-manifest-out ${sample}.plot_manifest.tsv
    """
}

process BUILD_COHORT_SUMMARY {
    tag "cohort_summary"
    publishDir "${params.outdir}/reports", mode: 'copy', pattern: 'cohort_*'
    publishDir "${params.outdir}/tables", mode: 'copy', pattern: 'cohort_*'

    input:
    path sample_summaries
    path gene_tables

    output:
    path 'cohort_sample_summary.tsv', emit: sample_summary
    path 'cohort_sample_summary.csv', emit: sample_summary_csv
    path 'cohort_gene_coverage_heatmap.html', emit: heatmap

    script:
    def sampleArgs = sample_summaries.collect { it.toString() }.join(' ')
    def geneArgs = gene_tables.collect { it.toString() }.join(' ')
    """
    ${params.python_executable} ${projectDir}/bin/build_cohort_summary.py \
        --sample-summaries ${sampleArgs} \
        --gene-tables ${geneArgs} \
        --summary-tsv cohort_sample_summary.tsv \
        --summary-csv cohort_sample_summary.csv \
        --heatmap-out cohort_gene_coverage_heatmap.html
    """
}
