# Dados de teste

Este diretório contém um conjunto sintético mínimo para exercício do pipeline.

## Como gerar os arquivos de entrada

```bash
cd test_data
bash generate_inputs.sh
```

Após a execução, será criada a pasta `test_data/input/` com:

- `sample1.bam`
- `sample1.bam.bai`
- `sample1.vcf.gz`
- `sample1.vcf.gz.tbi`
- `panel.bed`

## Observações

- Os dados são artificiais e servem apenas para validar fluxo, parsing e relatório.
- O VCF de teste contém uma variante PASS em região alvo e uma variante fora do painel.
