# QC logic

As regras de decisão são configuráveis via parâmetros Nextflow e serializadas em `results/qc/qc_rules.json`.

## Regras principais

- `pct_mapped >= min_mapped_reads_pct`
- `duplication_rate <= max_duplication_rate`
- `mean_target_coverage >= min_target_mean_coverage`
- `pct_target_ge_100x >= min_pct_target_100x`
- `pct_target_ge_200x >= min_pct_target_200x`
- `genes_with_inadequate_coverage <= max_critical_genes_below_threshold`
- `pass_variant_fraction >= min_variant_pass_fraction`
- `weak_variant_fraction <= max_weak_variant_fraction`

## Três classes de saída

- `LIBERAR`: sem flags `FAIL` ou `WARN`
- `LIBERAR COM RESSALVAS`: sem `FAIL`, mas com pelo menos uma `WARN`
- `NÃO LIBERAR`: uma ou mais regras em `FAIL`

## Tolerância para `WARN`

O pipeline usa `warning_tolerance_fraction` para definir uma zona intermediária:

- Para métricas com limite mínimo, valores entre `threshold * (1 - tolerance)` e `threshold` geram `WARN`
- Para métricas com limite máximo, valores entre `threshold` e `threshold * (1 + tolerance)` geram `WARN`
- Valores fora dessa faixa geram `FAIL`

## Definição de uniformidade

Uniformidade é definida como o percentual de bases alvo com profundidade maior ou igual a `0.2 * mean_target_coverage`.
