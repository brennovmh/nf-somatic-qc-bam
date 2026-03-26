#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$SCRIPT_DIR"

mkdir -p input
cp panel.bed input/panel.bed
cp sample1.vcf input/sample1.vcf

samtools faidx reference.fa
samtools view -bS sample1.sam | samtools sort -o input/sample1.bam
samtools index input/sample1.bam
bgzip -c input/sample1.vcf > input/sample1.vcf.gz
tabix -f -p vcf input/sample1.vcf.gz
rm -f input/sample1.vcf

echo 'Arquivos sintéticos gerados em test_data/input'
