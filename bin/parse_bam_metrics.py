#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path


def parse_flagstat(path: Path) -> dict:
    metrics = {}
    with open(path, 'r', encoding='utf-8') as handle:
        for line in handle:
            stripped = line.strip()
            if 'in total' in stripped:
                metrics['total_reads'] = int(stripped.split()[0])
            elif 'mapped (' in stripped and 'primary mapped' not in stripped:
                metrics['mapped_reads'] = int(stripped.split()[0])
                pct_match = re.search(r'\(([\d.]+)%', stripped)
                metrics['pct_mapped'] = float(pct_match.group(1)) if pct_match else None
            elif 'properly paired' in stripped:
                metrics['properly_paired_reads'] = int(stripped.split()[0])
                pct_match = re.search(r'\(([\d.]+)%', stripped)
                metrics['pct_properly_paired'] = float(pct_match.group(1)) if pct_match else None
    return metrics


def parse_stats(path: Path):
    metrics = {}
    mapq_rows = []
    total_reads = None
    with open(path, 'r', encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('SN\t'):
                _, label, value = line.rstrip('\n').split('\t')[:3]
                label = label.rstrip(':')
                if label == 'raw total sequences':
                    total_reads = int(float(value))
                elif label == 'reads duplicated':
                    metrics['duplicate_reads'] = int(float(value))
                elif label == 'average quality':
                    metrics['average_base_quality'] = float(value)
            elif line.startswith('MAPQ\t'):
                _, mapq, count = line.rstrip('\n').split('\t')
                count_int = int(count)
                mapq_rows.append({'mapq': int(mapq), 'count': count_int})
    for row in mapq_rows:
        row['pct_total'] = round((row['count'] / total_reads) * 100, 4) if total_reads else ''
    return metrics, mapq_rows


def write_table(path: Path, rows, fieldnames) -> None:
    with open(path, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--flagstat', required=True)
    parser.add_argument('--stats', required=True)
    parser.add_argument('--metrics-out', required=True)
    parser.add_argument('--mapq-out', required=True)
    args = parser.parse_args()

    flagstat = parse_flagstat(Path(args.flagstat))
    stats, mapq_rows = parse_stats(Path(args.stats))
    merged = {'sample_id': args.sample}
    merged.update(flagstat)
    merged.update(stats)
    if merged.get('duplicate_reads') is not None and merged.get('total_reads'):
        merged['duplication_rate'] = round(float(merged['duplicate_reads']) / float(merged['total_reads']), 6)

    write_table(Path(args.metrics_out), [merged], list(merged.keys()))
    mapq_out_rows = [{'sample_id': args.sample, **row} for row in mapq_rows]
    write_table(Path(args.mapq_out), mapq_out_rows, ['sample_id', 'mapq', 'count', 'pct_total'])


if __name__ == '__main__':
    main()
