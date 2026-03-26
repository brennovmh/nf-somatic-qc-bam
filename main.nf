nextflow.enable.dsl=2

include { VALIDATE_INPUTS }       from './modules/local/input_validation'
include { NORMALIZE_BED }         from './modules/local/panel'
include { WRITE_QC_RULES }        from './modules/local/common'
include { COLLECT_TOOL_VERSIONS } from './modules/local/common'
include { BUILD_COHORT_SUMMARY }  from './modules/local/reporting'
include { SAMPLE_QC }             from './subworkflows/local/sample_qc'


def findBaiFor(File bamFile) {
    def candidates = [ new File(bamFile.parentFile, bamFile.name + '.bai'), new File(bamFile.parentFile, bamFile.name.replaceFirst(/\.bam$/, '.bai')) ]
    def hit = candidates.find { it.exists() && it.length() > 0 }
    if( !hit ) error "Indice BAI nao encontrado para ${bamFile.name}"
    return hit
}

def inferVcfSample(File vcfFile) {
    if( vcfFile.name.endsWith('.vcf.gz') ) return vcfFile.name[0..-8]
    if( vcfFile.name.endsWith('.vcf') ) return vcfFile.name[0..-5]
    error "VCF com extensao nao suportada: ${vcfFile.name}"
}

def discoverSampleTuples(String inputDir, String samplePattern) {
    def input = new File(inputDir)
    if( !input.exists() || !input.isDirectory() ) error "--input_dir nao e um diretorio valido: ${inputDir}"

    def pattern = java.util.regex.Pattern.compile(samplePattern)
    def bamFiles = (input.listFiles() ?: [] as File[]).findAll { it.name.endsWith('.bam') && it.isFile() }.sort { it.name }
    def vcfFiles = (input.listFiles() ?: [] as File[]).findAll { (it.name.endsWith('.vcf') || it.name.endsWith('.vcf.gz')) && it.isFile() }.sort { it.name }

    if( !bamFiles ) error "Nenhum BAM encontrado em ${inputDir}"
    if( !vcfFiles ) error "Nenhum VCF encontrado em ${inputDir}"

    def bamMap = [:]
    bamFiles.each { bam ->
        def matcher = pattern.matcher(bam.name)
        if( !matcher.matches() ) error "BAM com nome fora do padrao esperado (${samplePattern}): ${bam.name}"
        def sample = matcher.group(1)
        bamMap[sample] = [bam: bam, bai: findBaiFor(bam)]
    }

    def vcfMap = [:]
    vcfFiles.each { vcf ->
        def sample = inferVcfSample(vcf)
        vcfMap[sample] = [vcf: vcf]
    }

    def missingVcf = (bamMap.keySet() - vcfMap.keySet()).sort()
    def missingBam = (vcfMap.keySet() - bamMap.keySet()).sort()
    if( missingVcf || missingBam ) {
        def parts = []
        if( missingVcf ) parts << "Sem VCF correspondente para: ${missingVcf.join(', ')}"
        if( missingBam ) parts << "Sem BAM correspondente para: ${missingBam.join(', ')}"
        error "Inconsistencia de amostras entre BAM e VCF. ${parts.join(' | ')}"
    }

    return bamMap.keySet().sort().collect { sample ->
        def bamInfo = bamMap[sample]
        def vcfInfo = vcfMap[sample]
        tuple(sample as String, file(bamInfo.bam), file(bamInfo.bai), file(vcfInfo.vcf), '')
    }
}


def qcRulesMap() {
    [
        min_mapped_reads_pct               : params.min_mapped_reads_pct,
        max_duplication_rate               : params.max_duplication_rate,
        min_target_mean_coverage           : params.min_target_mean_coverage,
        min_pct_target_100x                : params.min_pct_target_100x,
        min_pct_target_200x                : params.min_pct_target_200x,
        max_critical_genes_below_threshold : params.max_critical_genes_below_threshold,
        min_variant_pass_fraction          : params.min_variant_pass_fraction,
        min_variant_dp                     : params.min_variant_dp,
        min_variant_vaf                    : params.min_variant_vaf,
        max_weak_variant_fraction          : params.max_weak_variant_fraction,
        gene_primary_threshold             : params.gene_primary_threshold,
        min_gene_pct_primary               : params.min_gene_pct_primary,
        interval_primary_threshold         : params.interval_primary_threshold,
        min_interval_pct_primary           : params.min_interval_pct_primary,
        strand_bias_alert_threshold        : params.strand_bias_alert_threshold,
        warning_tolerance_fraction         : params.warning_tolerance_fraction
    ]
}

workflow {
    if( !params.input_dir ) {
        error "Parametro obrigatorio ausente: --input_dir"
    }

    def inputDirResolved = file(params.input_dir).toString()
    def bedResolved = params.bed ? file(params.bed).toString() : ''

    validated = VALIDATE_INPUTS(inputDirResolved, bedResolved, params.sample_pattern)
    normalized = NORMALIZE_BED(validated.panel_bed)
    rules_json = WRITE_QC_RULES(Channel.value(qcRulesMap()))
    COLLECT_TOOL_VERSIONS()

    samples_ch = Channel.fromList(discoverSampleTuples(inputDirResolved, params.sample_pattern as String))

    sample_qc = SAMPLE_QC(
        samples_ch,
        normalized.normalized_bed,
        normalized.panel_genes,
        normalized.panel_metadata,
        rules_json.rules_json
    )

    BUILD_COHORT_SUMMARY(
        sample_qc.sample_summary.collect(),
        sample_qc.gene_coverage.collect()
    )
}

workflow.onComplete {
    log.info "Pipeline finalizado. Status: ${workflow.success ? 'OK' : 'FALHA'}"
    log.info "Saidas em: ${params.outdir}"
}
