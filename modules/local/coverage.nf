process SAMTOOLS_TARGET_DEPTH {
    tag "${sample}"
    publishDir "${params.outdir}/coverage", mode: 'copy', pattern: '*.target.depth.tsv.gz'

    input:
    tuple val(sample), path(bam), path(bai), path(vcf), val(vcf_index)
    path bed_file

    output:
    tuple val(sample), path("${sample}.target.depth.tsv.gz"), emit: depth

    script:
    """
    samtools depth -@ ${task.cpus > 0 ? task.cpus - 1 : 0} -aa -d 0 -b ${bed_file} ${bam} | bgzip -c > ${sample}.target.depth.tsv.gz
    """
}

process AGGREGATE_COVERAGE {
    tag "${sample}"
    publishDir "${params.outdir}/coverage", mode: 'copy', pattern: '*.coverage_summary.tsv'
    publishDir "${params.outdir}/coverage", mode: 'copy', pattern: '*.gene_coverage.tsv'
    publishDir "${params.outdir}/coverage", mode: 'copy', pattern: '*.interval_coverage.tsv'
    publishDir "${params.outdir}/coverage", mode: 'copy', pattern: '*.low_coverage_*.tsv'

    input:
    tuple val(sample), path(depth_gz)
    path panel_genes

    output:
    tuple val(sample), path("${sample}.coverage_summary.tsv"), emit: coverage_summary
    tuple val(sample), path("${sample}.gene_coverage.tsv"), emit: gene_coverage
    tuple val(sample), path("${sample}.interval_coverage.tsv"), emit: interval_coverage
    tuple val(sample), path("${sample}.low_coverage_genes.tsv"), emit: low_coverage_genes
    tuple val(sample), path("${sample}.low_coverage_intervals.tsv"), emit: low_coverage_intervals

    script:
    """
    ${params.python_executable} ${projectDir}/bin/aggregate_coverage.py \
        --sample ${sample} \
        --depth ${depth_gz} \
        --panel-genes ${panel_genes} \
        --threshold-list 20,50,100,200,500 \
        --gene-primary-threshold ${params.gene_primary_threshold} \
        --min-gene-primary-pct ${params.min_gene_pct_primary} \
        --interval-primary-threshold ${params.interval_primary_threshold} \
        --min-interval-primary-pct ${params.min_interval_pct_primary} \
        --summary-out ${sample}.coverage_summary.tsv \
        --gene-out ${sample}.gene_coverage.tsv \
        --interval-out ${sample}.interval_coverage.tsv \
        --low-genes-out ${sample}.low_coverage_genes.tsv \
        --low-intervals-out ${sample}.low_coverage_intervals.tsv
    """
}
