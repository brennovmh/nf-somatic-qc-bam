process ESTIMATE_TUMOR_PURITY {
    tag "${sample}"
    publishDir "${params.outdir}/tables", mode: 'copy', pattern: '*.tumor_purity_estimate.tsv'

    input:
    tuple val(sample), path(variants)

    output:
    tuple val(sample), path("${sample}.tumor_purity_estimate.tsv"), emit: summary

    script:
    """
    ${params.python_executable} ${projectDir}/bin/estimate_tumor_purity.py \
        --sample ${sample} \
        --variants ${variants} \
        --min-dp ${params.tumor_purity_min_dp} \
        --min-vaf ${params.tumor_purity_min_vaf} \
        --max-vaf ${params.tumor_purity_max_vaf} \
        --min-variants ${params.tumor_purity_min_variants} \
        --bin-width ${params.tumor_purity_bin_width} \
        --enabled ${params.enable_tumor_purity} \
        --out ${sample}.tumor_purity_estimate.tsv
    """
}
