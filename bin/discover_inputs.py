#!/usr/bin/env python3
import argparse
import json
import re
import shutil
import sys
from pathlib import Path


def fail(message: str) -> None:
    print(f"ERRO: {message}", file=sys.stderr)
    sys.exit(1)


def non_empty(path: Path, label: str) -> None:
    if not path.exists():
        fail(f"{label} nao encontrado: {path}")
    if path.stat().st_size == 0:
        fail(f"{label} vazio: {path}")


def sample_from_bam(path: Path, pattern: re.Pattern) -> str:
    match = pattern.match(path.name)
    if not match:
        fail(f"BAM com nome fora do padrao esperado ({pattern.pattern}): {path.name}")
    sample = match.group(1)
    if not sample:
        fail(f"Nao foi possivel inferir o sample_id para {path.name}")
    return sample


def sample_from_vcf(path: Path) -> str:
    if path.name.endswith('.vcf.gz'):
        return path.name[:-7]
    if path.name.endswith('.vcf'):
        return path.name[:-4]
    fail(f"VCF com extensao nao suportada: {path.name}")


def find_bed(input_dir: Path, bed_override: str) -> Path:
    if bed_override:
        bed = Path(bed_override).resolve()
        non_empty(bed, 'BED informado em --bed')
        return bed
    beds = sorted(input_dir.glob('*.bed'))
    if not beds:
        fail('Nenhum arquivo BED encontrado no diretorio de entrada e --bed nao foi informado')
    if len(beds) > 1:
        fail('Mais de um arquivo BED encontrado. Informe explicitamente --bed')
    non_empty(beds[0], 'BED do painel')
    return beds[0].resolve()


def find_bai(bam: Path) -> Path:
    candidates = [Path(f'{bam}.bai'), bam.with_suffix('.bai')]
    for candidate in candidates:
        if candidate.exists():
            non_empty(candidate, 'Indice BAI')
            return candidate.resolve()
    fail(f'Indice BAI nao encontrado para {bam.name}')


def find_vcf_index(vcf: Path) -> str:
    if vcf.suffix == '.gz':
        for suffix in ('.tbi', '.csi'):
            candidate = Path(f'{vcf}{suffix}')
            if candidate.exists():
                non_empty(candidate, 'Indice VCF')
                return str(candidate.resolve())
    return ''


def main() -> None:
    parser = argparse.ArgumentParser(description='Descobre BAM/VCF e valida consistencia de entrada')
    parser.add_argument('--input-dir', required=True)
    parser.add_argument('--bed', default='')
    parser.add_argument('--sample-pattern', default=r'^(.*)\.bam$')
    parser.add_argument('--samplesheet', required=True)
    parser.add_argument('--panel-copy', required=True)
    parser.add_argument('--summary', required=True)
    args = parser.parse_args()

    input_dir = Path(args.input_dir).resolve()
    if not input_dir.is_dir():
        fail(f'--input_dir nao e um diretorio valido: {input_dir}')

    pattern = re.compile(args.sample_pattern)
    bam_files = sorted(p for p in input_dir.glob('*.bam') if p.is_file())
    if not bam_files:
        fail(f'Nenhum BAM encontrado em {input_dir}')

    vcf_files = sorted([p for p in input_dir.glob('*.vcf')] + [p for p in input_dir.glob('*.vcf.gz')])
    if not vcf_files:
        fail(f'Nenhum VCF encontrado em {input_dir}')

    bed = find_bed(input_dir, args.bed)

    bam_map = {}
    for bam in bam_files:
        non_empty(bam, 'BAM')
        sample = sample_from_bam(bam, pattern)
        if sample in bam_map:
            fail(f'Sample duplicado via BAM: {sample}')
        bam_map[sample] = {'bam': str(bam.resolve()), 'bai': str(find_bai(bam))}

    vcf_map = {}
    for vcf in vcf_files:
        non_empty(vcf, 'VCF')
        sample = sample_from_vcf(vcf)
        if sample in vcf_map:
            fail(f'Sample duplicado via VCF: {sample}')
        vcf_map[sample] = {'vcf': str(vcf.resolve()), 'vcf_index': find_vcf_index(vcf)}

    bam_samples = set(bam_map)
    vcf_samples = set(vcf_map)
    missing_vcf = sorted(bam_samples - vcf_samples)
    missing_bam = sorted(vcf_samples - bam_samples)
    if missing_vcf or missing_bam:
        messages = []
        if missing_vcf:
            messages.append(f"Sem VCF correspondente para: {', '.join(missing_vcf)}")
        if missing_bam:
            messages.append(f"Sem BAM correspondente para: {', '.join(missing_bam)}")
        fail('Inconsistencia de amostras entre BAM e VCF. ' + ' | '.join(messages))

    samples = []
    for sample in sorted(bam_samples):
        row = {'sample_id': sample}
        row.update(bam_map[sample])
        row.update(vcf_map[sample])
        samples.append(row)

    with open(args.samplesheet, 'w', encoding='utf-8') as handle:
        handle.write('sample_id\tbam\tbai\tvcf\tvcf_index\n')
        for row in samples:
            handle.write(f"{row['sample_id']}\t{row['bam']}\t{row['bai']}\t{row['vcf']}\t{row['vcf_index']}\n")

    shutil.copyfile(bed, args.panel_copy)

    summary = {'input_dir': str(input_dir), 'bed': str(bed), 'n_samples': len(samples), 'samples': [row['sample_id'] for row in samples]}
    with open(args.summary, 'w', encoding='utf-8') as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)


if __name__ == '__main__':
    main()
