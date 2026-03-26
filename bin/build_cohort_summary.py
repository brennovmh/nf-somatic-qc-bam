#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd
import plotly.express as px


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample-summaries', nargs='+', required=True)
    parser.add_argument('--gene-tables', nargs='+', required=True)
    parser.add_argument('--summary-tsv', required=True)
    parser.add_argument('--summary-csv', required=True)
    parser.add_argument('--heatmap-out', required=True)
    args = parser.parse_args()

    summaries = [pd.read_csv(path, sep='\t') for path in args.sample_summaries if Path(path).exists()]
    genes = [pd.read_csv(path, sep='\t') for path in args.gene_tables if Path(path).exists()]

    summary_df = pd.concat(summaries, ignore_index=True) if summaries else pd.DataFrame()
    summary_df.to_csv(args.summary_tsv, sep='\t', index=False)
    summary_df.to_csv(args.summary_csv, index=False)

    if genes:
        gene_df = pd.concat(genes, ignore_index=True)
        heatmap_df = gene_df.pivot_table(index='gene', columns='sample_id', values='mean_coverage', aggfunc='mean')
        fig = px.imshow(heatmap_df.fillna(0), aspect='auto', color_continuous_scale='Viridis', title='Heatmap de cobertura média por gene x amostra')
        Path(args.heatmap_out).write_text(fig.to_html(full_html=True, include_plotlyjs='cdn'), encoding='utf-8')
    else:
        Path(args.heatmap_out).write_text('<html><body><p>Sem dados para heatmap.</p></body></html>', encoding='utf-8')


if __name__ == '__main__':
    main()
