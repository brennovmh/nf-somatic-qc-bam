process PARSE_VCF {
    tag "${sample}"
    publishDir "${params.outdir}/variants", mode: 'copy', pattern: '*.vcf_summary.tsv'
    publishDir "${params.outdir}/variants", mode: 'copy', pattern: '*.variants*.tsv'
    publishDir "${params.outdir}/variants", mode: 'copy', pattern: '*.variant_gene_counts.tsv'

    input:
    tuple val(sample), path(bam), path(bai), path(vcf), val(vcf_index)
    path panel_genes

    output:
    tuple val(sample), path("${sample}.vcf_summary.tsv"), emit: vcf_summary
    tuple val(sample), path("${sample}.variants.tsv"), emit: variants
    tuple val(sample), path("${sample}.variants_pass.tsv"), emit: pass_variants
    tuple val(sample), path("${sample}.variants_nonpass.tsv"), emit: nonpass_variants
    tuple val(sample), path("${sample}.variant_gene_counts.tsv"), emit: variant_gene_counts

    script:
    """
    ${params.python_executable} ${projectDir}/bin/parse_vcf.py \
        --sample ${sample} \
        --vcf ${vcf} \
        --panel-genes ${panel_genes} \
        --min-variant-dp ${params.min_variant_dp} \
        --min-variant-vaf ${params.min_variant_vaf} \
        --strand-bias-threshold ${params.strand_bias_alert_threshold} \
        --summary-out ${sample}.vcf_summary.tsv \
        --variants-out ${sample}.variants.tsv \
        --pass-out ${sample}.variants_pass.tsv \
        --nonpass-out ${sample}.variants_nonpass.tsv \
        --gene-counts-out ${sample}.variant_gene_counts.tsv
    """
}
