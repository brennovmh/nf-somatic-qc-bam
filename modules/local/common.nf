process WRITE_QC_RULES {
    tag "qc_rules"
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: 'qc_rules.json'

    input:
    val rules

    output:
    path 'qc_rules.json', emit: rules_json

    script:
    def rulesJson = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(rules))
    """
    cat <<'JSON' > qc_rules.json
${rulesJson}
JSON
    """
}

process COLLECT_TOOL_VERSIONS {
    tag "tool_versions"
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: 'tool_versions.yml'

    output:
    path 'tool_versions.yml', emit: versions

    script:
    """
    {
      echo "python: \"\$(${params.python_executable} --version 2>&1 | awk '{print \$2}')\""
      echo "nextflow: \"${workflow.nextflow.version}\""
      echo "samtools: \"\$(samtools --version 2>/dev/null | head -n 1 | awk '{print \$2}' || echo not_available)\""
      echo "mosdepth: \"\$(mosdepth --version 2>/dev/null | awk '{print \$2}' || echo not_available)\""
      echo "bcftools: \"\$(bcftools --version 2>/dev/null | head -n 1 | awk '{print \$2}' || echo not_available)\""
      echo "bedtools: \"\$(bedtools --version 2>/dev/null | sed 's/bedtools v//' || echo not_available)\""
      echo "picard: \"\$(picard MarkDuplicates --version 2>/dev/null || echo not_available)\""
    } > tool_versions.yml
    """
}
