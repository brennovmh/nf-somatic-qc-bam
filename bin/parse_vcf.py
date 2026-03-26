#!/usr/bin/env python3
import argparse
import csv
import gzip
import math
import statistics
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path


def open_maybe_gzip(path):
    path = str(path)
    if path.endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8')
    return open(path, 'r', encoding='utf-8')


def load_panel(path):
    chrom_map = defaultdict(list)
    with open(path, 'r', encoding='utf-8') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            row['start'] = int(row['start'])
            row['end'] = int(row['end'])
            chrom_map[row['chrom']].append(row)
    starts = {}
    for chrom, rows in chrom_map.items():
        rows.sort(key=lambda x: x['start'])
        starts[chrom] = [row['start'] for row in rows]
    return chrom_map, starts


def overlaps(panel_rows, starts, chrom, pos0):
    chrom_rows = panel_rows.get(chrom, [])
    if not chrom_rows:
        return []
    idx = bisect_right(starts[chrom], pos0)
    hits = []
    for row in chrom_rows[max(0, idx - 20): idx + 20]:
        if row['start'] <= pos0 < row['end']:
            hits.append(row)
    return hits


def infer_variant_type(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return 'SNV'
    if len(ref) == len(alt):
        return 'MNV'
    return 'INDEL'


def mean_or_blank(values):
    values = [v for v in values if v is not None]
    return round(sum(values) / len(values), 6) if values else ''


def median_or_blank(values):
    values = [v for v in values if v is not None]
    return round(statistics.median(values), 6) if values else ''


def parse_info_field(info_str):
    info = {}
    if not info_str or info_str == '.':
        return info
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info[key] = value
        else:
            info[item] = True
    return info


def parse_sample_field(format_str, sample_str):
    if not format_str or not sample_str or format_str == '.' or sample_str == '.':
        return {}
    keys = format_str.split(':')
    values = sample_str.split(':')
    data = {}
    for idx, key in enumerate(keys):
        data[key] = values[idx] if idx < len(values) else ''
    return data


def parse_number(value, cast=float):
    if value in (None, '', '.'):
        return None
    try:
        return cast(value)
    except Exception:
        return None


def parse_list_numbers(value, cast=float):
    if value in (None, '', '.'):
        return []
    out = []
    for item in str(value).split(','):
        parsed = parse_number(item, cast=cast)
        out.append(parsed)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--vcf', required=True)
    parser.add_argument('--panel-genes', required=True)
    parser.add_argument('--min-variant-dp', type=float, required=True)
    parser.add_argument('--min-variant-vaf', type=float, required=True)
    parser.add_argument('--strand-bias-threshold', type=float, required=True)
    parser.add_argument('--summary-out', required=True)
    parser.add_argument('--variants-out', required=True)
    parser.add_argument('--pass-out', required=True)
    parser.add_argument('--nonpass-out', required=True)
    parser.add_argument('--gene-counts-out', required=True)
    args = parser.parse_args()

    panel_rows, starts = load_panel(Path(args.panel_genes))
    vcf_sample_name = None

    rows = []
    gene_counts = defaultdict(lambda: {'sample_id': args.sample, 'gene': None, 'total_variants': 0, 'pass_variants': 0})
    ti = 0
    tv = 0
    pass_dp = []
    pass_vaf = []
    low_dp_variants = 0
    low_vaf_variants = 0
    outside_panel = 0
    strand_bias_high = 0

    with open_maybe_gzip(args.vcf) as handle:
        for line in handle:
            if not line.strip():
                continue
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                header = line.rstrip('\n').split('\t')
                if len(header) > 9:
                    vcf_sample_name = header[9]
                continue

            fields = line.rstrip('\n').split('\t')
            if len(fields) < 8:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt_values = fields[4].split(',') if fields[4] not in ('', '.') else []
            qual = fields[5]
            filter_value = fields[6]
            info_map = parse_info_field(fields[7])
            format_str = fields[8] if len(fields) > 8 else ''
            sample_str = fields[9] if len(fields) > 9 else ''
            sample_data = parse_sample_field(format_str, sample_str)

            pos0 = pos - 1
            matches = overlaps(panel_rows, starts, chrom, pos0)
            genes = sorted({row['gene'] for row in matches if row['gene'] != 'NA'})
            interval_ids = sorted({row['interval_id'] for row in matches})

            dp_sample = parse_number(sample_data.get('DP'), cast=float)
            dp_info = parse_number(info_map.get('DP'), cast=float)
            af_sample = parse_list_numbers(sample_data.get('AF'), cast=float)
            af_info = parse_list_numbers(info_map.get('AF'), cast=float)
            ad_sample = parse_list_numbers(sample_data.get('AD'), cast=float)
            mq = parse_number(info_map.get('MQ'), cast=float)
            sb_candidates = parse_list_numbers(sample_data.get('SB'), cast=float)
            if not sb_candidates:
                sb_candidates = parse_list_numbers(info_map.get('SB'), cast=float)
            numeric_sb = None
            valid_sb = [x for x in sb_candidates if x is not None]
            if valid_sb:
                numeric_sb = max(valid_sb)

            for alt_index, alt in enumerate(alt_values):
                variant_type = infer_variant_type(ref, alt)
                dp = dp_sample if dp_sample is not None else dp_info
                vaf = None
                if alt_index < len(af_sample) and af_sample[alt_index] is not None:
                    vaf = af_sample[alt_index]
                elif alt_index < len(af_info) and af_info[alt_index] is not None:
                    vaf = af_info[alt_index]
                elif dp and len(ad_sample) > alt_index + 1 and ad_sample[alt_index + 1] is not None:
                    vaf = ad_sample[alt_index + 1] / dp if dp else None

                in_panel = bool(matches)
                filter_keys = [] if filter_value in ('', '.', 'PASS') else filter_value.split(';')
                is_pass = (not filter_keys) or filter_value == 'PASS'
                filters = 'PASS' if is_pass else ';'.join(filter_keys)

                if variant_type == 'SNV' and is_pass:
                    transition_pairs = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
                    if (ref.upper(), alt.upper()) in transition_pairs:
                        ti += 1
                    else:
                        tv += 1
                if is_pass and dp is not None:
                    pass_dp.append(float(dp))
                if is_pass and vaf is not None:
                    pass_vaf.append(float(vaf))
                if dp is not None and float(dp) < args.min_variant_dp:
                    low_dp_variants += 1
                if vaf is not None and float(vaf) < args.min_variant_vaf:
                    low_vaf_variants += 1
                if not in_panel:
                    outside_panel += 1
                if numeric_sb is not None and numeric_sb >= args.strand_bias_threshold:
                    strand_bias_high += 1

                gene_label = ','.join(genes) if genes else 'NA'
                interval_label = ','.join(interval_ids) if interval_ids else 'NA'
                row = {
                    'sample_id': args.sample,
                    'vcf_sample_name': vcf_sample_name or '',
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'variant_type': variant_type,
                    'qual': qual if qual != '.' else '',
                    'filter': filters,
                    'is_pass': is_pass,
                    'dp': dp if dp is not None else '',
                    'vaf': round(float(vaf), 6) if vaf is not None and not math.isnan(float(vaf)) else '',
                    'mq': mq if mq is not None else '',
                    'strand_bias': numeric_sb if numeric_sb is not None else '',
                    'gene': gene_label,
                    'interval_id': interval_label,
                    'in_panel': in_panel,
                    'weak_support': (
                        (dp is not None and float(dp) < args.min_variant_dp)
                        or (vaf is not None and float(vaf) < args.min_variant_vaf)
                        or (numeric_sb is not None and numeric_sb >= args.strand_bias_threshold)
                        or (not in_panel)
                    ),
                }
                rows.append(row)
                for gene in genes or ['NA']:
                    gene_counts[gene]['gene'] = gene
                    gene_counts[gene]['total_variants'] += 1
                    if is_pass:
                        gene_counts[gene]['pass_variants'] += 1

    fieldnames = list(rows[0].keys()) if rows else [
        'sample_id', 'vcf_sample_name', 'chrom', 'pos', 'ref', 'alt', 'variant_type', 'qual', 'filter',
        'is_pass', 'dp', 'vaf', 'mq', 'strand_bias', 'gene', 'interval_id', 'in_panel', 'weak_support'
    ]
    with open(args.variants_out, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    pass_rows = [row for row in rows if row['is_pass']]
    nonpass_rows = [row for row in rows if not row['is_pass']]
    for output_path, subset in [(args.pass_out, pass_rows), (args.nonpass_out, nonpass_rows)]:
        with open(output_path, 'w', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, delimiter='\t', fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(subset)

    gene_count_rows = sorted(gene_counts.values(), key=lambda row: row['gene'] or 'ZZZ')
    with open(args.gene_counts_out, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=['sample_id', 'gene', 'total_variants', 'pass_variants'])
        writer.writeheader()
        writer.writerows(gene_count_rows)

    total = len(rows)
    passed = len(pass_rows)
    summary = {
        'sample_id': args.sample,
        'vcf_sample_name': vcf_sample_name or '',
        'total_variants': total,
        'pass_variants': passed,
        'filtered_variants': len(nonpass_rows),
        'pct_pass': round((passed / total) * 100, 4) if total else 0.0,
        'snv_count': sum(1 for row in rows if row['variant_type'] == 'SNV'),
        'indel_count': sum(1 for row in rows if row['variant_type'] == 'INDEL'),
        'mnv_count': sum(1 for row in rows if row['variant_type'] == 'MNV'),
        'mean_pass_dp': mean_or_blank(pass_dp),
        'median_pass_dp': median_or_blank(pass_dp),
        'mean_pass_vaf': mean_or_blank(pass_vaf),
        'median_pass_vaf': median_or_blank(pass_vaf),
        'titv_ratio': round(ti / tv, 4) if tv else '',
        'low_dp_variants': low_dp_variants,
        'low_vaf_variants': low_vaf_variants,
        'variants_outside_panel': outside_panel,
        'strand_bias_high_variants': strand_bias_high,
    }
    with open(args.summary_out, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=list(summary.keys()))
        writer.writeheader()
        writer.writerow(summary)


if __name__ == '__main__':
    main()
