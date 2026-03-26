process NORMALIZE_BED {
    tag "normalize_panel"
    publishDir "${params.outdir}/tables", mode: 'copy', pattern: 'panel*'

    input:
    path bed_file

    output:
    path 'panel.normalized.bed', emit: normalized_bed
    path 'panel_genes.tsv', emit: panel_genes
    path 'panel_metadata.json', emit: panel_metadata

    script:
    """
    ${params.python_executable} ${projectDir}/bin/normalize_bed.py \
        --bed ${bed_file} \
        --normalized-bed panel.normalized.bed \
        --panel-genes panel_genes.tsv \
        --metadata panel_metadata.json
    """
}
