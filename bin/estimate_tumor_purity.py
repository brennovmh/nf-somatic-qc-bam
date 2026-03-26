#!/usr/bin/env python3
import argparse
import csv
import math
import statistics


def parse_float(value):
    if value in (None, '', '.'):
        return None
    try:
        return float(value)
    except Exception:
        return None


def truthy(value):
    return str(value).strip().lower() in {'true', '1', 'yes'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--variants', required=True)
    parser.add_argument('--min-dp', type=float, required=True)
    parser.add_argument('--min-vaf', type=float, required=True)
    parser.add_argument('--max-vaf', type=float, required=True)
    parser.add_argument('--min-variants', type=int, required=True)
    parser.add_argument('--bin-width', type=float, required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--enabled', required=True)
    args = parser.parse_args()

    enabled = str(args.enabled).strip().lower() in {'true', '1', 'yes'}

    eligible_vafs = []
    pass_snv_count = 0
    total_variants = 0
    total_pass = 0

    with open(args.variants, 'r', encoding='utf-8') as handle:
        reader = csv.DictReader(handle, delimiter='	')
        for row in reader:
            total_variants += 1
            is_pass = truthy(row.get('is_pass'))
            variant_type = str(row.get('variant_type', '')).upper()
            in_panel = truthy(row.get('in_panel'))
            weak_support = truthy(row.get('weak_support'))
            dp = parse_float(row.get('dp'))
            vaf = parse_float(row.get('vaf'))
            if is_pass:
                total_pass += 1
            if is_pass and variant_type == 'SNV':
                pass_snv_count += 1
            if not is_pass or variant_type != 'SNV' or not in_panel or weak_support:
                continue
            if dp is None or dp < args.min_dp:
                continue
            if vaf is None or vaf < args.min_vaf or vaf > args.max_vaf:
                continue
            eligible_vafs.append(vaf)

    method = 'Heuristica de pico dominante de VAF em SNVs PASS no painel'
    assumptions = 'Assume variantes clonais heterozigoticas em regioes diploides; pureza aproximada = 2 x VAF dominante'
    limitations = 'Estimativa exploratoria sem normal pareado; CNA, LOH, aneuploidia e subclonalidade podem enviesar o resultado'

    summary = {
        'sample_id': args.sample,
        'method': method,
        'enabled': enabled,
        'total_variants': total_variants,
        'pass_variants': total_pass,
        'pass_snv_variants': pass_snv_count,
        'eligible_variants': len(eligible_vafs),
        'min_dp': args.min_dp,
        'min_vaf': args.min_vaf,
        'max_vaf': args.max_vaf,
        'bin_width': args.bin_width,
        'dominant_vaf_peak': '',
        'peak_bin_start': '',
        'peak_bin_end': '',
        'peak_supporting_variants': 0,
        'peak_support_fraction': '',
        'estimated_tumor_purity': '',
        'confidence_class': 'NOT_ESTIMABLE',
        'status': 'DISABLED' if not enabled else 'NOT_ESTIMABLE',
        'summary_text': 'Modulo de pureza tumoral desabilitado.' if not enabled else 'Estimativa nao calculavel: variantes elegiveis insuficientes.',
        'assumptions': assumptions,
        'limitations': limitations,
    }

    if enabled and eligible_vafs and args.bin_width > 0:
        n_bins = max(1, int(math.ceil((args.max_vaf - args.min_vaf) / args.bin_width)))
        bins = [[] for _ in range(n_bins)]
        for value in eligible_vafs:
            idx = int((value - args.min_vaf) / args.bin_width)
            if idx >= n_bins:
                idx = n_bins - 1
            bins[idx].append(value)
        peak_idx, peak_values = max(enumerate(bins), key=lambda item: (len(item[1]), -item[0]))
        peak_start = args.min_vaf + (peak_idx * args.bin_width)
        peak_end = min(args.max_vaf, peak_start + args.bin_width)
        peak_vaf = statistics.median(peak_values) if peak_values else None
        peak_fraction = (len(peak_values) / len(eligible_vafs)) if eligible_vafs else 0.0
        estimated = min(1.0, max(0.0, 2 * peak_vaf)) if peak_vaf is not None else None

        summary['dominant_vaf_peak'] = round(peak_vaf, 6) if peak_vaf is not None else ''
        summary['peak_bin_start'] = round(peak_start, 6)
        summary['peak_bin_end'] = round(peak_end, 6)
        summary['peak_supporting_variants'] = len(peak_values)
        summary['peak_support_fraction'] = round(peak_fraction, 6)
        summary['estimated_tumor_purity'] = round(estimated, 6) if estimated is not None else ''

        if len(eligible_vafs) < args.min_variants:
            summary['confidence_class'] = 'NOT_ESTIMABLE'
            summary['status'] = 'NOT_ESTIMABLE'
            summary['summary_text'] = 'Estimativa nao calculavel: variantes elegiveis insuficientes.'
        elif len(eligible_vafs) < max(5, args.min_variants + 1) or peak_fraction < 0.35:
            summary['confidence_class'] = 'LOW_CONFIDENCE'
            summary['status'] = 'ESTIMATED'
            summary['summary_text'] = (
                f'Estimativa heuristica com baixa confianca: pureza aproximada {summary["estimated_tumor_purity"]}, '
                f'baseada em {len(eligible_vafs)} SNVs elegiveis e pico dominante de VAF em {summary["dominant_vaf_peak"]}.'
            )
        else:
            summary['confidence_class'] = 'MODERATE_CONFIDENCE'
            summary['status'] = 'ESTIMATED'
            summary['summary_text'] = (
                f'Estimativa heuristica de pureza aproximada {summary["estimated_tumor_purity"]}, '
                f'baseada em {len(eligible_vafs)} SNVs elegiveis e pico dominante de VAF em {summary["dominant_vaf_peak"]}.'
            )

    fieldnames = [
        'sample_id', 'method', 'enabled', 'total_variants', 'pass_variants', 'pass_snv_variants',
        'eligible_variants', 'min_dp', 'min_vaf', 'max_vaf', 'bin_width', 'dominant_vaf_peak',
        'peak_bin_start', 'peak_bin_end', 'peak_supporting_variants', 'peak_support_fraction',
        'estimated_tumor_purity', 'confidence_class', 'status', 'summary_text', 'assumptions', 'limitations'
    ]
    with open(args.out, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, delimiter='	', fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(summary)


if __name__ == '__main__':
    main()
