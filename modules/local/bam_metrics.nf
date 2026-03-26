process SAMTOOLS_FLAGSTAT {
    tag "${sample}"
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*.flagstat.txt'

    input:
    tuple val(sample), path(bam), path(bai), path(vcf), val(vcf_index)

    output:
    tuple val(sample), path("${sample}.flagstat.txt"), emit: txt

    script:
    """
    samtools flagstat -@ ${task.cpus > 0 ? task.cpus - 1 : 0} ${bam} > ${sample}.flagstat.txt
    """
}

process SAMTOOLS_STATS {
    tag "${sample}"
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*.samtools.stats.txt'

    input:
    tuple val(sample), path(bam), path(bai), path(vcf), val(vcf_index)

    output:
    tuple val(sample), path("${sample}.samtools.stats.txt"), emit: txt

    script:
    """
    samtools stats -@ ${task.cpus > 0 ? task.cpus - 1 : 0} ${bam} > ${sample}.samtools.stats.txt
    """
}

process PARSE_BAM_METRICS {
    tag "${sample}"
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*.bam_metrics.tsv'
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*.mapq_distribution.tsv'

    input:
    tuple val(sample), path(flagstat_txt), path(stats_txt)

    output:
    tuple val(sample), path("${sample}.bam_metrics.tsv"), emit: bam_metrics
    tuple val(sample), path("${sample}.mapq_distribution.tsv"), emit: mapq_distribution

    script:
    """
    ${params.python_executable} ${projectDir}/bin/parse_bam_metrics.py \
        --sample ${sample} \
        --flagstat ${flagstat_txt} \
        --stats ${stats_txt} \
        --metrics-out ${sample}.bam_metrics.tsv \
        --mapq-out ${sample}.mapq_distribution.tsv
    """
}
