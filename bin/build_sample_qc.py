#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import pandas as pd


def warn_fail_for_min(observed, threshold, tolerance):
    if pd.isna(observed):
        return 'WARN'
    if observed >= threshold:
        return 'PASS'
    if observed >= threshold * (1 - tolerance):
        return 'WARN'
    return 'FAIL'


def warn_fail_for_max(observed, threshold, tolerance):
    if pd.isna(observed):
        return 'WARN'
    if observed <= threshold:
        return 'PASS'
    if observed <= threshold * (1 + tolerance):
        return 'WARN'
    return 'FAIL'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--bam-metrics', required=True)
    parser.add_argument('--coverage-summary', required=True)
    parser.add_argument('--gene-coverage', required=True)
    parser.add_argument('--interval-coverage', required=True)
    parser.add_argument('--vcf-summary', required=True)
    parser.add_argument('--variants', required=True)
    parser.add_argument('--tumor-purity', required=True)
    parser.add_argument('--panel-genes', required=True)
    parser.add_argument('--rules', required=True)
    parser.add_argument('--summary-tsv', required=True)
    parser.add_argument('--summary-csv', required=True)
    parser.add_argument('--flags-out', required=True)
    parser.add_argument('--classification-out', required=True)
    parser.add_argument('--executive-out', required=True)
    parser.add_argument('--panel-status-out', required=True)
    parser.add_argument('--tumor-purity-out', required=True)
    args = parser.parse_args()

    rules = json.loads(Path(args.rules).read_text(encoding='utf-8'))
    bam = pd.read_csv(args.bam_metrics, sep='\t')
    coverage = pd.read_csv(args.coverage_summary, sep='\t')
    genes = pd.read_csv(args.gene_coverage, sep='\t') if Path(args.gene_coverage).stat().st_size else pd.DataFrame()
    intervals = pd.read_csv(args.interval_coverage, sep='\t') if Path(args.interval_coverage).stat().st_size else pd.DataFrame()
    vcf = pd.read_csv(args.vcf_summary, sep='\t')
    variants = pd.read_csv(args.variants, sep='\t') if Path(args.variants).stat().st_size else pd.DataFrame()
    tumor_purity = pd.read_csv(args.tumor_purity, sep='\t') if Path(args.tumor_purity).stat().st_size else pd.DataFrame()
    panel = pd.read_csv(args.panel_genes, sep='\t')

    purity_for_row = tumor_purity.drop(columns=['sample_id'], errors='ignore') if not tumor_purity.empty else pd.DataFrame()
    row = pd.concat([bam, coverage, vcf.drop(columns=['sample_id'], errors='ignore'), purity_for_row], axis=1)
    row = row.loc[:, ~row.columns.duplicated()].copy()
    tolerance = float(rules['warning_tolerance_fraction'])

    low_genes = genes[genes['status'] != 'ADEQUATE'] if not genes.empty and 'status' in genes.columns else pd.DataFrame()
    low_intervals = intervals[intervals['status'] != 'ADEQUATE'] if not intervals.empty and 'status' in intervals.columns else pd.DataFrame()

    weak_variants = variants[
        (variants['weak_support'] == True) | (variants['weak_support'].astype(str).str.lower() == 'true')
    ] if not variants.empty and 'weak_support' in variants.columns else pd.DataFrame()
    weak_variant_fraction = (len(weak_variants) / len(variants)) if len(variants) else 0.0

    panel_status = []
    if not genes.empty and 'gene' in panel.columns:
        pass_variants = variants[variants['is_pass'].astype(str).str.lower() == 'true'].copy() if not variants.empty else pd.DataFrame()
        gene_variant_counts = (
            pass_variants.assign(gene=pass_variants['gene'].fillna('NA').astype(str).str.split(','))
            .explode('gene')
            .groupby('gene')
            .size()
            .to_dict()
        ) if not pass_variants.empty else {}
        for gene, gene_rows in panel[panel['gene'] != 'NA'].groupby('gene'):
            gene_cov = genes[genes['gene'] == gene]
            coverage_adequate = bool(not gene_cov.empty and gene_cov.iloc[0]['status'] == 'ADEQUATE')
            pass_count = int(gene_variant_counts.get(gene, 0))
            panel_status.append({
                'sample_id': args.sample,
                'gene': gene,
                'num_intervals': int(len(gene_rows)),
                'target_size': int(gene_rows['target_size'].sum()),
                'coverage_adequate': coverage_adequate,
                'gene_status': gene_cov.iloc[0]['status'] if not gene_cov.empty else 'NA',
                'pass_variant_count': pass_count,
                'has_pass_variant': pass_count > 0,
                'adequate_but_no_variant': coverage_adequate and pass_count == 0,
            })

    flags = []

    def add_flag(metric, observed, comparator, threshold, status, message):
        flags.append({
            'sample_id': args.sample,
            'metric': metric,
            'observed': observed,
            'comparator': comparator,
            'threshold': threshold,
            'status': status,
            'message': message,
        })

    pct_mapped = float(row.at[0, 'pct_mapped'])
    add_flag('pct_mapped', pct_mapped, '>=', rules['min_mapped_reads_pct'], warn_fail_for_min(pct_mapped, float(rules['min_mapped_reads_pct']), tolerance), 'Percentual de reads mapeados')

    dup_rate = float(row.at[0, 'duplication_rate']) if 'duplication_rate' in row.columns else float('nan')
    add_flag('duplication_rate', dup_rate, '<=', rules['max_duplication_rate'], warn_fail_for_max(dup_rate, float(rules['max_duplication_rate']), tolerance), 'Taxa de duplicacao')

    mean_cov = float(row.at[0, 'mean_target_coverage'])
    add_flag('mean_target_coverage', mean_cov, '>=', rules['min_target_mean_coverage'], warn_fail_for_min(mean_cov, float(rules['min_target_mean_coverage']), tolerance), 'Cobertura media alvo')

    pct_100 = float(row.at[0, 'pct_target_ge_100x'])
    add_flag('pct_target_ge_100x', pct_100, '>=', rules['min_pct_target_100x'], warn_fail_for_min(pct_100, float(rules['min_pct_target_100x']), tolerance), 'Percentual de bases alvo >=100x')

    pct_200 = float(row.at[0, 'pct_target_ge_200x'])
    add_flag('pct_target_ge_200x', pct_200, '>=', rules['min_pct_target_200x'], warn_fail_for_min(pct_200, float(rules['min_pct_target_200x']), tolerance), 'Percentual de bases alvo >=200x')

    low_genes_count = int(len(low_genes))
    add_flag('critical_genes_low_coverage', low_genes_count, '<=', rules['max_critical_genes_below_threshold'], warn_fail_for_max(low_genes_count, int(rules['max_critical_genes_below_threshold']), tolerance), 'Genes criticos com cobertura inadequada')

    pass_fraction = float(row.at[0, 'pct_pass']) / 100.0 if 'pct_pass' in row.columns else 0.0
    add_flag('pass_variant_fraction', pass_fraction, '>=', rules['min_variant_pass_fraction'], warn_fail_for_min(pass_fraction, float(rules['min_variant_pass_fraction']), tolerance), 'Fracao de variantes PASS')

    add_flag('weak_variant_fraction', round(weak_variant_fraction, 6), '<=', rules['max_weak_variant_fraction'], warn_fail_for_max(weak_variant_fraction, float(rules['max_weak_variant_fraction']), tolerance), 'Fracao de variantes com suporte tecnico fraco')

    fail_count = sum(1 for flag in flags if flag['status'] == 'FAIL')
    warn_count = sum(1 for flag in flags if flag['status'] == 'WARN')
    if fail_count > 0:
        final_classification = 'NÃO LIBERAR'
        final_statement = 'Dados não aptos para análise terciária'
    elif warn_count > 0:
        final_classification = 'LIBERAR COM RESSALVAS'
        final_statement = 'Dados aptos com ressalvas'
    else:
        final_classification = 'LIBERAR'
        final_statement = 'Dados aptos para análise terciária'

    failed_messages = [flag['message'] for flag in flags if flag['status'] == 'FAIL']
    warning_messages = [flag['message'] for flag in flags if flag['status'] == 'WARN']
    rationale_parts = []
    if failed_messages:
        rationale_parts.append('Falhas críticas: ' + '; '.join(failed_messages))
    if warning_messages:
        rationale_parts.append('Alertas: ' + '; '.join(warning_messages))
    if not rationale_parts:
        rationale_parts.append('Todas as métricas principais atenderam aos critérios configurados.')
    rationale = ' '.join(rationale_parts)

    row['final_classification'] = final_classification
    row['final_statement'] = final_statement
    row['fail_count'] = fail_count
    row['warn_count'] = warn_count
    row['genes_with_inadequate_coverage'] = low_genes_count
    row['intervals_with_inadequate_coverage'] = len(low_intervals)
    row['weak_variant_fraction'] = round(weak_variant_fraction, 6)
    row['rationale'] = rationale

    row.to_csv(args.summary_tsv, sep='\t', index=False)
    row.to_csv(args.summary_csv, index=False)

    pd.DataFrame(flags).to_csv(args.flags_out, sep='\t', index=False)
    pd.DataFrame([{
        'sample_id': args.sample,
        'final_classification': final_classification,
        'final_statement': final_statement,
        'rationale': rationale,
    }]).to_csv(args.classification_out, sep='\t', index=False)

    with open(args.executive_out, 'w', encoding='utf-8') as handle:
        handle.write(f'# Resumo executivo - {args.sample}\n\n')
        handle.write(f'- Classificação final: **{final_classification}**\n')
        handle.write(f'- Conclusão: **{final_statement}**\n')
        handle.write(f'- Cobertura média alvo: **{mean_cov:.2f}x**\n')
        handle.write(f'- Percentual de bases alvo >=100x: **{pct_100:.2f}%**\n')
        handle.write(f'- Genes com cobertura inadequada: **{low_genes_count}**\n')
        handle.write(f'- Variantes PASS: **{int(row.at[0, "pass_variants"]) if "pass_variants" in row.columns else 0}**\n')
        if not tumor_purity.empty:
            purity_row = tumor_purity.iloc[0]
            purity_value = purity_row['estimated_tumor_purity']
            purity_text = purity_value if pd.notna(purity_value) and purity_value != '' else 'NA'
            handle.write(f'- Pureza tumoral heurística: **{purity_text}**\n')
            handle.write(f'- Confiança da pureza: **{purity_row["confidence_class"]}**\n')
        handle.write(f'- Justificativa: {rationale}\n')

    panel_status_df = pd.DataFrame(panel_status, columns=['sample_id', 'gene', 'num_intervals', 'target_size', 'coverage_adequate', 'gene_status', 'pass_variant_count', 'has_pass_variant', 'adequate_but_no_variant'])
    panel_status_df.to_csv(args.panel_status_out, sep='\t', index=False)
    tumor_purity.to_csv(args.tumor_purity_out, sep='\t', index=False)


if __name__ == '__main__':
    main()
