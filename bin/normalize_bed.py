#!/usr/bin/env python3
import argparse
import csv
import json
import re
from pathlib import Path


def clean_gene(value: str) -> str:
    text = (value or '').strip()
    if not text:
        return 'NA'
    text = re.split(r'[|,;]', text)[0].strip()
    return text or 'NA'


def main() -> None:
    parser = argparse.ArgumentParser(description='Normaliza BED do painel')
    parser.add_argument('--bed', required=True)
    parser.add_argument('--normalized-bed', required=True)
    parser.add_argument('--panel-genes', required=True)
    parser.add_argument('--metadata', required=True)
    args = parser.parse_args()

    bed_path = Path(args.bed)
    if not bed_path.exists() or bed_path.stat().st_size == 0:
        raise SystemExit(f'BED invalido: {bed_path}')

    intervals = []
    gene_detected = False
    with open(bed_path, 'r', encoding='utf-8') as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith(('#', 'track', 'browser')):
                continue
            fields = raw.rstrip('\n').split('\t')
            if len(fields) < 3:
                raise SystemExit(f'Linha BED invalida: {raw.rstrip()}')
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            if end <= start:
                raise SystemExit(f'Intervalo BED invalido: {raw.rstrip()}')
            gene = clean_gene(fields[3]) if len(fields) >= 4 else 'NA'
            if gene != 'NA':
                gene_detected = True
            intervals.append((chrom, start, end, gene))

    if not intervals:
        raise SystemExit('BED sem intervalos validos')

    intervals.sort(key=lambda x: (x[0], x[1], x[2], x[3]))

    with open(args.normalized_bed, 'w', encoding='utf-8') as bed_out, open(args.panel_genes, 'w', encoding='utf-8', newline='') as table_out:
        writer = csv.writer(table_out, delimiter='\t')
        writer.writerow(['chrom', 'start', 'end', 'gene', 'interval_id', 'label', 'target_size'])
        for index, (chrom, start, end, gene) in enumerate(intervals, start=1):
            interval_id = f'INT{index:06d}'
            label = gene if gene != 'NA' else interval_id
            bed_out.write(f'{chrom}\t{start}\t{end}\t{label}\t{interval_id}\n')
            writer.writerow([chrom, start, end, gene, interval_id, label, end - start])

    metadata = {
        'source_bed': str(bed_path.resolve()),
        'interval_count': len(intervals),
        'has_gene_mapping': gene_detected,
        'gene_inference_rule': 'BED coluna 4, com primeira etiqueta antes de | , ou ;' if gene_detected else 'Sem gene na coluna 4; agregacao por gene indisponivel'
    }
    with open(args.metadata, 'w', encoding='utf-8') as handle:
        json.dump(metadata, handle, indent=2, ensure_ascii=False)


if __name__ == '__main__':
    main()
