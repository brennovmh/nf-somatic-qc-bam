#!/usr/bin/env python3
import argparse
import csv
import gzip
import statistics
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path


def read_panel(path: Path):
    panel_rows = []
    with open(path, 'r', encoding='utf-8') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            row['start'] = int(row['start'])
            row['end'] = int(row['end'])
            row['target_size'] = int(row['target_size'])
            panel_rows.append(row)
    panel_rows.sort(key=lambda r: (r['chrom'], r['start'], r['end']))
    return panel_rows


def build_chrom_index(panel_rows):
    chrom_index = defaultdict(list)
    for row in panel_rows:
        chrom_index[row['chrom']].append(row)
    starts = {chrom: [item['start'] for item in items] for chrom, items in chrom_index.items()}
    return chrom_index, starts


def overlapping_intervals(chrom_rows, chrom_starts, pos0):
    idx = bisect_right(chrom_starts, pos0)
    hits = []
    for candidate in chrom_rows[max(0, idx - 20): idx + 20]:
        if candidate['start'] <= pos0 < candidate['end']:
            hits.append(candidate['interval_id'])
    return hits


def summarize_depths(depths, thresholds):
    if not depths:
        return {
            'mean_coverage': 0.0,
            'median_coverage': 0.0,
            'min_coverage': 0,
            'max_coverage': 0,
            **{f'pct_ge_{t}x': 0.0 for t in thresholds},
        }
    summary = {
        'mean_coverage': round(sum(depths) / len(depths), 4),
        'median_coverage': round(statistics.median(depths), 4),
        'min_coverage': min(depths),
        'max_coverage': max(depths),
    }
    for threshold in thresholds:
        count = sum(1 for depth in depths if depth >= threshold)
        summary[f'pct_ge_{threshold}x'] = round((count / len(depths)) * 100, 4)
    return summary


def write_tsv(path: Path, rows, fieldnames):
    with open(path, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--depth', required=True)
    parser.add_argument('--panel-genes', required=True)
    parser.add_argument('--threshold-list', required=True)
    parser.add_argument('--gene-primary-threshold', type=int, required=True)
    parser.add_argument('--min-gene-primary-pct', type=float, required=True)
    parser.add_argument('--interval-primary-threshold', type=int, required=True)
    parser.add_argument('--min-interval-primary-pct', type=float, required=True)
    parser.add_argument('--summary-out', required=True)
    parser.add_argument('--gene-out', required=True)
    parser.add_argument('--interval-out', required=True)
    parser.add_argument('--low-genes-out', required=True)
    parser.add_argument('--low-intervals-out', required=True)
    args = parser.parse_args()

    thresholds = [int(x) for x in args.threshold_list.split(',')]
    panel_rows = read_panel(Path(args.panel_genes))
    chrom_index, start_index = build_chrom_index(panel_rows)
    interval_by_id = {row['interval_id']: row for row in panel_rows}
    interval_depths = {row['interval_id']: [] for row in panel_rows}
    all_depths = []

    with gzip.open(args.depth, 'rt', encoding='utf-8') as handle:
        for line in handle:
            chrom, pos, depth = line.rstrip('\n').split('\t')
            pos0 = int(pos) - 1
            depth_value = int(depth)
            all_depths.append(depth_value)
            chrom_rows = chrom_index.get(chrom)
            if not chrom_rows:
                continue
            hits = overlapping_intervals(chrom_rows, start_index[chrom], pos0)
            for interval_id in hits:
                interval_depths[interval_id].append(depth_value)

    interval_rows = []
    low_interval_rows = []
    gene_groups = defaultdict(list)
    for interval_id, row in interval_by_id.items():
        depths = interval_depths[interval_id]
        stats = summarize_depths(depths, thresholds)
        interval_row = {
            'sample_id': args.sample,
            'gene': row['gene'],
            'interval_id': interval_id,
            'chrom': row['chrom'],
            'start': row['start'],
            'end': row['end'],
            'target_size': row['target_size'],
            'mean_coverage': stats['mean_coverage'],
            'median_coverage': stats['median_coverage'],
            'min_coverage': stats['min_coverage'],
            'max_coverage': stats['max_coverage'],
        }
        for threshold in thresholds:
            interval_row[f'pct_ge_{threshold}x'] = stats[f'pct_ge_{threshold}x']
        interval_row['status'] = 'LOW_COVERAGE' if interval_row[f'pct_ge_{args.interval_primary_threshold}x'] < args.min_interval_primary_pct else 'ADEQUATE'
        interval_rows.append(interval_row)
        gene_groups[row['gene']].append((interval_row, depths))
        if interval_row['status'] != 'ADEQUATE':
            low_interval_rows.append(interval_row)

    gene_rows = []
    low_gene_rows = []
    for gene, entries in sorted(gene_groups.items()):
        if gene == 'NA':
            continue
        combined_depths = []
        total_size = 0
        for interval_row, depths in entries:
            combined_depths.extend(depths)
            total_size += interval_row['target_size']
        stats = summarize_depths(combined_depths, thresholds)
        gene_row = {
            'sample_id': args.sample,
            'gene': gene,
            'num_intervals': len(entries),
            'target_size': total_size,
            'mean_coverage': stats['mean_coverage'],
            'median_coverage': stats['median_coverage'],
            'min_coverage': stats['min_coverage'],
            'max_coverage': stats['max_coverage'],
        }
        for threshold in thresholds:
            gene_row[f'pct_ge_{threshold}x'] = stats[f'pct_ge_{threshold}x']
        gene_row['status'] = 'INADEQUATE' if gene_row[f'pct_ge_{args.gene_primary_threshold}x'] < args.min_gene_primary_pct else 'ADEQUATE'
        gene_rows.append(gene_row)
        if gene_row['status'] != 'ADEQUATE':
            low_gene_rows.append(gene_row)

    global_stats = summarize_depths(all_depths, thresholds)
    mean_cov = global_stats['mean_coverage']
    uniform_count = sum(1 for depth in all_depths if depth >= 0.2 * mean_cov) if all_depths else 0
    coverage_summary = {
        'sample_id': args.sample,
        'target_bases': len(all_depths),
        'mean_target_coverage': global_stats['mean_coverage'],
        'median_target_coverage': global_stats['median_coverage'],
        'min_target_coverage': global_stats['min_coverage'],
        'max_target_coverage': global_stats['max_coverage'],
        'uniformity_pct_ge_0.2x_mean': round((uniform_count / len(all_depths)) * 100, 4) if all_depths else 0.0,
        'genes_with_low_coverage': len(low_gene_rows),
        'intervals_with_low_coverage': len(low_interval_rows),
    }
    for threshold in thresholds:
        coverage_summary[f'pct_target_ge_{threshold}x'] = global_stats[f'pct_ge_{threshold}x']

    interval_fields = list(interval_rows[0].keys()) if interval_rows else [
        'sample_id', 'gene', 'interval_id', 'chrom', 'start', 'end', 'target_size', 'mean_coverage',
        'median_coverage', 'min_coverage', 'max_coverage', 'pct_ge_20x', 'pct_ge_50x', 'pct_ge_100x',
        'pct_ge_200x', 'pct_ge_500x', 'status'
    ]
    gene_fields = list(gene_rows[0].keys()) if gene_rows else [
        'sample_id', 'gene', 'num_intervals', 'target_size', 'mean_coverage', 'median_coverage',
        'min_coverage', 'max_coverage', 'pct_ge_20x', 'pct_ge_50x', 'pct_ge_100x', 'pct_ge_200x',
        'pct_ge_500x', 'status'
    ]

    write_tsv(Path(args.summary_out), [coverage_summary], list(coverage_summary.keys()))
    write_tsv(Path(args.interval_out), interval_rows, interval_fields)
    write_tsv(Path(args.gene_out), gene_rows, gene_fields)
    write_tsv(Path(args.low_genes_out), low_gene_rows, gene_fields)
    write_tsv(Path(args.low_intervals_out), low_interval_rows, interval_fields)


if __name__ == '__main__':
    main()
