process VALIDATE_INPUTS {
    tag "discover_inputs"
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: 'discovery_*'

    input:
    val input_dir
    val bed_override
    val sample_pattern

    output:
    path 'samplesheet.tsv', emit: samplesheet
    path 'panel.bed', emit: panel_bed
    path 'discovery_summary.json', emit: discovery_summary

    script:
    """
    ${params.python_executable} ${projectDir}/bin/discover_inputs.py \
        --input-dir "${input_dir}" \
        --bed "${bed_override}" \
        --sample-pattern '${sample_pattern}' \
        --samplesheet samplesheet.tsv \
        --panel-copy panel.bed \
        --summary discovery_summary.json
    """
}
