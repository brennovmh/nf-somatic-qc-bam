#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import pandas as pd
import plotly.express as px


def df_to_html(df: pd.DataFrame, max_rows: int = 20) -> str:
    if df.empty:
        return '<p>Sem dados disponíveis.</p>'
    if max_rows is None or max_rows <= 0:
        return df.to_html(index=False, classes='table table-sm', border=0)
    return df.head(max_rows).to_html(index=False, classes='table table-sm', border=0)


def searchable_table_html(df: pd.DataFrame, table_id: str, placeholder: str) -> str:
    if df.empty:
        return '<p>Sem dados disponíveis.</p>'
    table_html = df.to_html(index=False, classes='table table-sm searchable-table', border=0, table_id=table_id)
    return f"""
<div class='search-box'>
  <input type='search' id='{table_id}_search' placeholder='{placeholder}' onkeyup="filterTable('{table_id}')" />
</div>
<div class='table-scroll'>
  {table_html}
</div>
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--sample-summary', required=True)
    parser.add_argument('--classification', required=True)
    parser.add_argument('--qc-flags', required=True)
    parser.add_argument('--gene-coverage', required=True)
    parser.add_argument('--interval-coverage', required=True)
    parser.add_argument('--variants', required=True)
    parser.add_argument('--mapq-distribution', required=True)
    parser.add_argument('--tumor-purity', required=True)
    parser.add_argument('--panel-variant-status', required=True)
    parser.add_argument('--panel-metadata', required=True)
    parser.add_argument('--html-out', required=True)
    parser.add_argument('--markdown-out', required=True)
    parser.add_argument('--plot-manifest-out', required=True)
    args = parser.parse_args()

    sample_summary = pd.read_csv(args.sample_summary, sep='\t')
    classification = pd.read_csv(args.classification, sep='\t')
    qc_flags = pd.read_csv(args.qc_flags, sep='\t')
    gene_coverage = pd.read_csv(args.gene_coverage, sep='\t') if Path(args.gene_coverage).stat().st_size else pd.DataFrame()
    interval_coverage = pd.read_csv(args.interval_coverage, sep='\t') if Path(args.interval_coverage).stat().st_size else pd.DataFrame()
    variants = pd.read_csv(args.variants, sep='\t') if Path(args.variants).stat().st_size else pd.DataFrame()
    mapq = pd.read_csv(args.mapq_distribution, sep='\t') if Path(args.mapq_distribution).stat().st_size else pd.DataFrame()
    tumor_purity = pd.read_csv(args.tumor_purity, sep='\t') if Path(args.tumor_purity).stat().st_size else pd.DataFrame()
    panel_status = pd.read_csv(args.panel_variant_status, sep='\t') if Path(args.panel_variant_status).stat().st_size else pd.DataFrame()
    panel_metadata = json.loads(Path(args.panel_metadata).read_text(encoding='utf-8'))

    summary_row = sample_summary.iloc[0]
    class_row = classification.iloc[0]

    badge_class = {
        'LIBERAR': 'badge-pass',
        'LIBERAR COM RESSALVAS': 'badge-warn',
        'NÃO LIBERAR': 'badge-fail',
    }.get(class_row['final_classification'], 'badge-warn')

    plot_manifest = []
    plots_html = []

    if not gene_coverage.empty:
        fig = px.bar(gene_coverage.sort_values('mean_coverage', ascending=False), x='gene', y='mean_coverage', color='status', title='Cobertura média por gene', labels={'mean_coverage': 'Cobertura média (x)', 'gene': 'Gene'})
        fig.update_xaxes(showticklabels=False, title_text='Genes do painel')
        plots_html.append(fig.to_html(full_html=False, include_plotlyjs='cdn'))
        plot_manifest.append({'sample_id': args.sample, 'plot': 'gene_mean_coverage', 'type': 'plotly'})

    if not interval_coverage.empty:
        fig = px.histogram(interval_coverage, x='mean_coverage', nbins=40, title='Distribuição de cobertura média por intervalo', labels={'mean_coverage': 'Cobertura média por intervalo (x)'})
        plots_html.append(fig.to_html(full_html=False, include_plotlyjs=False))
        plot_manifest.append({'sample_id': args.sample, 'plot': 'interval_coverage_distribution', 'type': 'plotly'})

    if not variants.empty and 'vaf' in variants.columns:
        variant_vaf = variants[variants['vaf'].astype(str) != ''].copy()
        if not variant_vaf.empty:
            variant_vaf['vaf'] = pd.to_numeric(variant_vaf['vaf'])
            fig = px.histogram(variant_vaf, x='vaf', nbins=30, color='filter', title='Distribuição de VAF', labels={'vaf': 'VAF'})
            plots_html.append(fig.to_html(full_html=False, include_plotlyjs=False))
            plot_manifest.append({'sample_id': args.sample, 'plot': 'vaf_distribution', 'type': 'plotly'})

    if not variants.empty and 'dp' in variants.columns:
        variant_dp = variants[variants['dp'].astype(str) != ''].copy()
        if not variant_dp.empty:
            variant_dp['dp'] = pd.to_numeric(variant_dp['dp'])
            fig = px.histogram(variant_dp, x='dp', nbins=30, color='filter', title='Distribuição de profundidade nas variantes', labels={'dp': 'DP'})
            plots_html.append(fig.to_html(full_html=False, include_plotlyjs=False))
            plot_manifest.append({'sample_id': args.sample, 'plot': 'variant_dp_distribution', 'type': 'plotly'})

    if not mapq.empty:
        fig = px.bar(mapq.sort_values('mapq'), x='mapq', y='count', title='Distribuição de MAPQ', labels={'mapq': 'MAPQ', 'count': 'Reads'})
        plots_html.append(fig.to_html(full_html=False, include_plotlyjs=False))
        plot_manifest.append({'sample_id': args.sample, 'plot': 'mapq_distribution', 'type': 'plotly'})

    failed_flags = qc_flags[qc_flags['status'] == 'FAIL']
    warned_flags = qc_flags[qc_flags['status'] == 'WARN']
    low_genes = gene_coverage[gene_coverage['status'] != 'ADEQUATE'] if not gene_coverage.empty else pd.DataFrame()
    purity_row = tumor_purity.iloc[0] if not tumor_purity.empty else None
    weak_variants = variants[(variants['weak_support'].astype(str).str.lower() == 'true')] if not variants.empty and 'weak_support' in variants.columns else pd.DataFrame()

    html = f"""<!DOCTYPE html>
<html lang='pt-BR'>
<head>
  <meta charset='utf-8' />
  <title>QC secundário - {args.sample}</title>
  <style>
    body {{ font-family: Helvetica, Arial, sans-serif; margin: 2rem auto; max-width: 1280px; color: #17202a; line-height: 1.5; }}
    h1, h2, h3 {{ color: #0f3b52; }}
    .hero {{ padding: 1.5rem; border-radius: 16px; background: linear-gradient(135deg, #eff7f6, #d9edf7); margin-bottom: 1.5rem; }}
    .badge {{ display: inline-block; padding: 0.35rem 0.7rem; border-radius: 999px; color: #fff; font-weight: 700; }}
    .badge-pass {{ background: #1f7a4d; }}
    .badge-warn {{ background: #c97c00; }}
    .badge-fail {{ background: #b42318; }}
    .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 1rem; margin: 1rem 0 1.5rem 0; }}
    .card {{ border: 1px solid #d0d7de; border-radius: 14px; padding: 1rem; background: #fff; }}
    .table {{ width: 100%; border-collapse: collapse; }}
    .table th, .table td {{ border-bottom: 1px solid #e5e7eb; padding: 0.45rem; text-align: left; font-size: 0.92rem; }}
    .table-scroll {{ max-height: 620px; overflow: auto; border: 1px solid #d0d7de; border-radius: 14px; background: #fff; }}
    .search-box {{ margin: 0.5rem 0 0.85rem 0; }}
    .search-box input {{ width: min(420px, 100%); padding: 0.7rem 0.9rem; border: 1px solid #cbd5e1; border-radius: 10px; font-size: 0.95rem; }}
    .muted {{ color: #4b5563; }}
  </style>
  <script>
    function filterTable(tableId) {{
      const input = document.getElementById(tableId + '_search');
      const filter = ((input && input.value) || '').toLowerCase();
      const table = document.getElementById(tableId);
      if (!table) return;
      const rows = table.getElementsByTagName('tr');
      for (let i = 1; i < rows.length; i++) {{
        const rowText = (rows[i].innerText || rows[i].textContent || '').toLowerCase();
        rows[i].style.display = rowText.indexOf(filter) > -1 ? '' : 'none';
      }}
    }}
  </script>
</head>
<body>
  <div class='hero'>
    <h1>Relatório técnico de QC secundário</h1>
    <p><strong>Amostra:</strong> {args.sample}</p>
    <p><span class='badge {badge_class}'>{class_row['final_classification']}</span></p>
    <p><strong>Conclusão:</strong> {class_row['final_statement']}</p>
    <p><strong>Justificativa:</strong> {class_row['rationale']}</p>
  </div>

  <h2>Resumo executivo</h2>
  <div class='grid'>
    <div class='card'><strong>Cobertura média alvo</strong><br />{summary_row['mean_target_coverage']:.2f}x</div>
    <div class='card'><strong>% alvo &gt;=100x</strong><br />{summary_row['pct_target_ge_100x']:.2f}%</div>
    <div class='card'><strong>% reads mapeados</strong><br />{summary_row['pct_mapped']:.2f}%</div>
    <div class='card'><strong>Duplicação</strong><br />{summary_row['duplication_rate']:.4f}</div>
    <div class='card'><strong>Variantes PASS</strong><br />{int(summary_row['pass_variants'])}</div>
    <div class='card'><strong>Genes inadequados</strong><br />{int(summary_row['genes_with_inadequate_coverage'])}</div>
  </div>

  <h2>Alertas técnicos</h2>
  <h3>Falhas críticas</h3>
  {df_to_html(failed_flags)}
  <h3>Ressalvas</h3>
  {df_to_html(warned_flags)}

  <h2>Métricas principais do BAM e cobertura</h2>
  {sample_summary.to_html(index=False, border=0, classes='table')}

  <h2>Estimativa heurística de pureza tumoral</h2>
  {df_to_html(tumor_purity, max_rows=10) if not tumor_purity.empty else '<p>Módulo de pureza tumoral desabilitado ou sem dados.</p>'}
  <p class='muted'>
    {purity_row['summary_text'] if purity_row is not None else 'Estimativa não disponível.'}
  </p>

  <h2>Métricas por gene</h2>
  {searchable_table_html(gene_coverage.sort_values(['status', 'gene'], ascending=[False, True]), 'gene_metrics_table', 'Pesquisar gene, status ou cobertura...') if not gene_coverage.empty else '<p>Sem agregação por gene.</p>'}

  <h2>Genes com cobertura inadequada</h2>
  {df_to_html(low_genes, max_rows=50)}

  <h2>Gráficos</h2>
  {''.join(plots_html)}

  <h2>Metadados do painel</h2>
  <pre>{json.dumps(panel_metadata, indent=2, ensure_ascii=False)}</pre>

  <p class='muted'>Este relatório consolida BAM, VCF e BED para apoiar a decisão automatizada de liberação para análise terciária.</p>
</body>
</html>
"""
    Path(args.html_out).write_text(html, encoding='utf-8')

    markdown = f"""# Relatório técnico de QC secundário - {args.sample}

## Status final

- Classificação: **{class_row['final_classification']}**
- Conclusão: **{class_row['final_statement']}**
- Justificativa: {class_row['rationale']}

## Métricas-chave

- Cobertura média alvo: **{summary_row['mean_target_coverage']:.2f}x**
- Percentual de bases alvo >=100x: **{summary_row['pct_target_ge_100x']:.2f}%**
- Percentual de reads mapeados: **{summary_row['pct_mapped']:.2f}%**
- Taxa de duplicação: **{summary_row['duplication_rate']:.4f}**
- Genes com cobertura inadequada: **{int(summary_row['genes_with_inadequate_coverage'])}**
- Pureza tumoral heurística: **{summary_row['estimated_tumor_purity'] if 'estimated_tumor_purity' in summary_row.index and pd.notna(summary_row['estimated_tumor_purity']) and summary_row['estimated_tumor_purity'] != '' else 'NA'}**
- Confiança da pureza: **{summary_row['confidence_class'] if 'confidence_class' in summary_row.index else 'NA'}**

## Observações

- Total de flags FAIL: **{len(failed_flags)}**
- Total de flags WARN: **{len(warned_flags)}**
- Intervalos com baixa cobertura: **{int(summary_row['intervals_with_inadequate_coverage'])}**
"""
    Path(args.markdown_out).write_text(markdown, encoding='utf-8')
    pd.DataFrame(plot_manifest).to_csv(args.plot_manifest_out, sep='\t', index=False)


if __name__ == '__main__':
    main()
