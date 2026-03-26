# qc-secundario-somático

Pipeline em Nextflow DSL2 para controle de qualidade técnico da análise secundária de NGS somático em painel alvo, integrando BAM, VCF e BED para classificar cada amostra como `LIBERAR`, `LIBERAR COM RESSALVAS` ou `NÃO LIBERAR`.

## Objetivo

Receber uma pasta contendo:

- `*.bam`
- `*.bam.bai` ou `*.bai`
- `*.vcf` ou `*.vcf.gz`
- índices `*.tbi` ou `*.csi` para VCF compactado, quando existirem
- `*.bed` do painel

E gerar:

- métricas globais de alinhamento
- métricas de cobertura no alvo
- métricas por intervalo e por gene
- métricas técnicas do VCF
- flags automatizadas de QC
- decisão final sobre aptidão para análise terciária
- relatório HTML robusto para uso laboratorial

## Estrutura do projeto

```text
qc-secundario/
├── main.nf
├── nextflow.config
├── params.example.json
├── README.md
├── assets/
│   └── README.md
├── bin/
│   ├── aggregate_coverage.py
│   ├── build_cohort_summary.py
│   ├── build_sample_qc.py
│   ├── discover_inputs.py
│   ├── normalize_bed.py
│   ├── parse_bam_metrics.py
│   ├── parse_vcf.py
│   └── render_report.py
├── conf/
│   ├── base.config
│   ├── conda.config
│   ├── docker.config
│   └── standard.config
├── docker/
│   └── Dockerfile
├── docs/
│   ├── dag.txt
│   └── qc_logic.md
├── envs/
│   └── qc-secundario.yml
├── modules/
│   └── local/
│       ├── bam_metrics.nf
│       ├── common.nf
│       ├── coverage.nf
│       ├── input_validation.nf
│       ├── panel.nf
│       ├── reporting.nf
│       └── vcf_metrics.nf
├── subworkflows/
│   └── local/
│       └── sample_qc.nf
└── test_data/
    ├── README.md
    ├── generate_inputs.sh
    ├── panel.bed
    ├── reference.fa
    ├── sample1.sam
    └── sample1.vcf
```

## Requisitos

### Perfil `standard`

Espera ferramentas disponíveis no ambiente local:

- `python3`
- `samtools`
- `bcftools`
- `bedtools`
- `bgzip`
- `tabix`

### Perfil `docker`

Use a imagem definida em `params.container_image`.

Build sugerido:

```bash
docker build -t qc-secundario:1.0.0 -f docker/Dockerfile .
```

## Execução

### Exemplo mínimo

```bash
nextflow run main.nf \
  -profile standard \
  --input_dir /caminho/para/input \
  --outdir results
```

### Exemplo com BED explícito

```bash
nextflow run main.nf \
  -profile standard \
  --input_dir /caminho/para/input \
  --bed /caminho/para/painel.bed \
  --outdir results
```

### Exemplo com Docker

```bash
nextflow run main.nf \
  -profile docker \
  --input_dir /caminho/para/input \
  --outdir results \
  --container_image qc-secundario:1.0.0
```

### Exemplo com parâmetros em JSON

```bash
nextflow run main.nf \
  -params-file params.example.json \
  -profile standard
```

## Estrutura de input esperada

```text
input/
├── sample1.bam
├── sample1.bam.bai
├── sample1.vcf.gz
├── sample1.vcf.gz.tbi
├── sample2.bam
├── sample2.bam.bai
├── sample2.vcf.gz
├── sample2.vcf.gz.tbi
└── painel.bed
```

## Detecção de amostras

- O pipeline infere o `sample_id` a partir do nome do BAM com `sample_pattern`
- O VCF correspondente deve ter o mesmo prefixo do BAM
- Exemplos:
  - `sample1.bam` combina com `sample1.vcf.gz`
  - `sample2.bam` combina com `sample2.vcf`
- Se houver inconsistência entre BAM e VCF, o pipeline falha com mensagem clara

## Ferramentas usadas

- `samtools flagstat`: contagens globais de reads e mapeamento
- `samtools stats`: estatísticas adicionais e distribuição de MAPQ
- `samtools depth -aa -b`: profundidade base a base nas regiões alvo
- `pysam`: parsing tolerante de VCF e extração dinâmica de INFO/FORMAT
- `pandas` e `plotly`: integração das tabelas, decisão QC e relatório HTML

## Racional técnico

### Por que o profile `standard` usa apenas `samtools depth`

- No modo `standard`, toda a cobertura é calculada a partir de `samtools depth`, reduzindo fragilidade ambiental
- Isso mantém mediana, min, max, uniformidade, percentuais por threshold e agregação por gene/intervalo sem depender de `mosdepth`

### Por que o parsing de VCF foi feito em Python em vez de somente `bcftools query`

`bcftools` é excelente para estatísticas globais, mas em painéis clínicos reais há heterogeneidade de campos `INFO/FORMAT`. O parser em Python foi mantido tolerante a ausência de campos opcionais e, no profile `standard`, evita dependências frágeis adicionais no ambiente.

### MultiQC

`MultiQC` não foi incluído como etapa principal porque o objetivo aqui é um relatório decisório integrado BAM + VCF + BED com regras laboratoriais. Se desejado, pode ser adicionado depois como saída complementar.

## Métricas calculadas

### BAM

- número total de reads
- reads mapeados e percentual mapeado
- reads propriamente pareados, quando aplicável
- taxa de duplicação estimada a partir de `samtools stats`
- distribuição de MAPQ
- cobertura média e mediana nas regiões alvo
- percentual de bases alvo >=20x, >=50x, >=100x, >=200x e >=500x
- uniformidade de cobertura
- cobertura por intervalo
- cobertura agregada por gene
- lista de genes e intervalos com cobertura insuficiente
- depth mean, median, min e max por gene

### VCF

- número total de variantes
- variantes PASS e filtradas
- contagem de SNV, INDEL e MNV
- distribuição de VAF
- distribuição de DP em variantes
- média e mediana de DP nas variantes PASS
- média e mediana de VAF nas variantes PASS
- razão Ti/Tv para SNVs PASS quando aplicável
- variantes com baixa profundidade
- variantes com VAF abaixo do mínimo configurado
- variantes com possível strand bias alto
- checagem se a variante está ou não dentro do painel

### Integração

- variantes fora do painel são sinalizadas
- genes do painel com cobertura inadequada são destacados
- genes com cobertura adequada, porém sem variantes PASS, são listados
- regiões tecnicamente inadequadas são reportadas para interpretação cuidadosa

## Glossário das métricas e conceitos-chave

A tabela abaixo resume os principais campos usados nas saídas tabulares e no relatório HTML. Ela serve como referência rápida para interpretação técnica e laboratorial.

| Campo / conceito | Origem principal | O que significa | Interpretação prática |
| --- | --- | --- | --- |
| `total_reads` | BAM / `samtools flagstat` | Total de reads observados no BAM | Volume bruto de sequenciamento/alinhamento processado |
| `mapped_reads` | BAM / `samtools flagstat` | Reads marcados como mapeados | Base para avaliar sucesso do alinhamento |
| `pct_mapped` | BAM / `samtools flagstat` | Percentual de reads mapeados em relação ao total | Valor baixo sugere problema de alinhamento, contaminação ou material ruim |
| `properly_paired_reads` | BAM / `samtools flagstat` | Reads em pares corretos, quando a biblioteca é paired-end | Queda relevante pode sugerir artefatos estruturais ou alinhamento ruim |
| `duplication_rate` | BAM / `samtools stats` | Fração estimada de reads duplicados | Valor alto reduz complexidade útil e pode impactar sensibilidade |
| `mean_target_coverage` | BAM + BED / `samtools depth` | Cobertura média nas bases alvo do painel | Métrica global de profundidade útil para decisão de liberação |
| `median_target_coverage` | BAM + BED / `samtools depth` | Cobertura mediana nas bases alvo | Menos sensível a regiões hiperamplificadas que a média |
| `pct_target_ge_20x`, `pct_target_ge_50x`, `pct_target_ge_100x`, `pct_target_ge_200x`, `pct_target_ge_500x` | BAM + BED / `samtools depth` | Percentual das bases alvo acima de cada threshold | Mede completude técnica do painel em thresholds clinicamente relevantes |
| `coverage_uniformity` | BAM + BED / `samtools depth` | Percentual de bases alvo com depth >= 0.2 x cobertura média alvo | Resume homogeneidade da cobertura; quanto maior, melhor |
| `genes_with_inadequate_coverage` | Cobertura agregada por gene | Número de genes abaixo do critério configurado | Ajuda a estimar impacto técnico na interpretação terciária |
| `intervals_with_inadequate_coverage` | Cobertura por intervalo | Número de intervalos BED abaixo do critério configurado | Mostra falhas localizadas que podem afetar hotspots/regiões críticas |
| `gene status` | Cobertura agregada por gene | Classificação do gene como `ADEQUATE` ou `INADEQUATE` | Resume se o gene está tecnicamente interpretável conforme a regra QC |
| `total_variants` | VCF | Total de variantes lidas do VCF | Volume bruto de chamadas antes de priorização por filtro |
| `pass_variants` | VCF | Variantes com `FILTER=PASS` | Conjunto principal de chamadas tecnicamente aceitáveis |
| `filtered_variants` | VCF | Variantes com qualquer filtro diferente de `PASS` | Chamadas que exigem cautela ou exclusão, dependendo do pipeline primário |
| `pass_fraction` | VCF | Proporção de variantes PASS sobre o total | Queda importante pode sugerir ruído técnico ou chamada instável |
| `snv_count` | VCF | Número de SNVs | Perfil geral do tipo de variante chamada |
| `indel_count` | VCF | Número de INDELs | Útil para ver se há excesso anômalo de indels ou ausência inesperada |
| `mnv_count` | VCF | Número de MNVs quando identificáveis | Campo auxiliar; pode ser zero dependendo do caller |
| `mean_pass_dp` / `median_pass_dp` | VCF | Média / mediana da profundidade nas variantes PASS | Mede suporte de profundidade das variantes liberáveis |
| `mean_pass_vaf` / `median_pass_vaf` | VCF | Média / mediana de VAF nas variantes PASS | Resume o comportamento alélico das chamadas aprovadas |
| `vaf` | VCF `FORMAT/INFO` | Fração alélica variante da chamada | Muito baixo pode indicar limite técnico, subclonalidade ou ruído |
| `dp` | VCF `FORMAT/INFO` | Profundidade associada à variante | Base para confiança mínima na chamada |
| `weak_support` | Integração VCF + regras QC | Flag para variante com profundidade, VAF ou suporte técnico fraco | Sinaliza variantes que não devem ser interpretadas sem cautela |
| `outside_panel` | VCF + BED | Indica se a variante está fora das regiões alvo do painel | Pode apontar ruído, padding inesperado ou inconsistência de painel |
| `strand_bias_high` | VCF `INFO/FORMAT` | Flag para viés de fita acima do limite configurado, quando o campo existe | Sinal de possível artefato técnico |
| `panel genes without PASS variants` | Integração BAM + VCF + BED | Genes com cobertura adequada, mas sem variantes PASS | Não é falha técnica por si só; serve para contextualização do painel |
| `final_classification` | Integração de regras QC | Resultado final: `LIBERAR`, `LIBERAR COM RESSALVAS` ou `NÃO LIBERAR` | Síntese automatizada da aptidão para seguir à análise terciária |
| `final_statement` | Integração de regras QC | Texto conclusivo do relatório | Frase final de uso operacional no laboratório |
| `rationale` | Integração de regras QC | Justificativa objetiva baseada nas flags disparadas | Explica por que a amostra foi ou não liberada |

## QC logic

Resumo em `docs/qc_logic.md`.

Regras configuráveis por parâmetro:

- `min_mapped_reads_pct`
- `max_duplication_rate`
- `min_target_mean_coverage`
- `min_pct_target_100x`
- `min_pct_target_200x`
- `max_critical_genes_below_threshold`
- `min_variant_pass_fraction`
- `min_variant_dp`
- `min_variant_vaf`
- `max_weak_variant_fraction`

Classificação:

- `LIBERAR`: sem flags `FAIL` ou `WARN`
- `LIBERAR COM RESSALVAS`: pelo menos uma `WARN`, mas sem `FAIL`
- `NÃO LIBERAR`: qualquer flag `FAIL`

## Definição de métricas

### Uniformidade

Uniformidade = percentual de bases alvo com profundidade maior ou igual a `0.2 * mean_target_coverage`.

### Genes insuficientemente cobertos

Por padrão, um gene é marcado como inadequado quando `% de bases do gene >= gene_primary_threshold` fica abaixo de `min_gene_pct_primary`.

### Intervalos críticos

Por padrão, um intervalo é marcado como inadequado quando `% de bases do intervalo >= interval_primary_threshold` fica abaixo de `min_interval_pct_primary`.

## Mapeamento BED para gene

- O pipeline tenta inferir gene a partir da coluna 4 do BED
- Quando a coluna 4 contém texto composto, usa a primeira etiqueta antes de `|`, `,` ou `;`
- Se não houver gene confiável, o pipeline continua em modo por intervalo
- A suposição é registrada em `panel_metadata.json`

## Saídas

```text
results/
├── coverage/
├── logs/
├── qc/
├── reports/
├── tables/
└── variants/
```

Arquivos principais por amostra:

- `tables/<sample>.sample_summary.tsv`
- `tables/<sample>.sample_summary.csv`
- `coverage/<sample>.gene_coverage.tsv`
- `coverage/<sample>.interval_coverage.tsv`
- `variants/<sample>.variants.tsv`
- `variants/<sample>.variants_pass.tsv`
- `variants/<sample>.variants_nonpass.tsv`
- `qc/<sample>.qc_flags.tsv`
- `qc/<sample>.classification.tsv`
- `tables/<sample>.tumor_purity_estimate.tsv`
- `reports/<sample>.qc_report.html`
- `reports/<sample>.qc_report.md`

Arquivos de coorte:

- `tables/cohort_sample_summary.tsv`
- `tables/cohort_sample_summary.csv`
- `reports/cohort_gene_coverage_heatmap.html`

### Como interpretar cada arquivo de saída

#### Arquivos por amostra

| Arquivo | Onde fica | O que contém | Como interpretar / quando usar |
| --- | --- | --- | --- |
| `sample_summary.tsv` / `sample_summary.csv` | `results/tables/` | Resumo consolidado da amostra com métricas BAM, cobertura, VCF e decisão integrada | É a melhor tabela para visão geral rápida e comparação entre amostras |
| `coverage_summary.tsv` | `results/coverage/` | Métricas globais de cobertura do painel, como média, mediana, uniformidade e percentuais por threshold | Use para avaliar se a cobertura global atende aos critérios mínimos do laboratório |
| `gene_coverage.tsv` | `results/coverage/` | Métricas agregadas por gene, incluindo cobertura média, mediana, min, max e percentuais por threshold | Arquivo principal para revisar genes críticos do painel e falhas localizadas por gene |
| `interval_coverage.tsv` | `results/coverage/` | Métricas por intervalo do BED | Use quando for necessário investigar exon, hotspot ou região alvo específica |
| `low_coverage_genes.tsv` | `results/coverage/` | Subconjunto dos genes classificados como inadequados | Ideal para revisão rápida dos genes que podem comprometer a análise terciária |
| `low_coverage_intervals.tsv` | `results/coverage/` | Subconjunto dos intervalos classificados como inadequados | Útil para identificar falhas técnicas regionais e regiões difíceis recorrentes |
| `variants.tsv` | `results/variants/` | Tabela completa de variantes extraídas do VCF com campos técnicos relevantes | É a visão mais abrangente das variantes para auditoria e investigação técnica |
| `variants_pass.tsv` | `results/variants/` | Apenas variantes com `FILTER=PASS` | Arquivo mais útil para revisar o conjunto principal de variantes aproveitáveis |
| `variants_nonpass.tsv` | `results/variants/` | Variantes filtradas ou não PASS | Use para investigar ruído, filtros excessivos ou chamadas marginais |
| `variant_gene_counts.tsv` | `results/variants/` | Contagem de variantes por gene, quando o mapeamento é possível | Ajuda a contextualizar genes do painel com e sem eventos detectados |
| `bam_metrics.tsv` | `results/qc/` | Métricas extraídas de `samtools flagstat` e `samtools stats` | Melhor arquivo para revisar alinhamento, mapeamento, pareamento e duplicação |
| `mapq_distribution.tsv` | `results/qc/` | Distribuição da qualidade de mapeamento | Use quando houver suspeita de alinhamento subótimo ou excesso de leituras de baixa confiança |
| `vcf_summary.tsv` | `results/variants/` | Resumo global do VCF com contagens, fração PASS, DP/VAF e Ti/Tv quando aplicável | Arquivo-chave para entender a qualidade geral das chamadas de variantes |
| `qc_flags.tsv` | `results/qc/` | Lista de flags `PASS`, `WARN` e `FAIL` disparadas pelas regras configuradas | Deve ser revisado sempre que a amostra vier com ressalva ou não liberação |
| `classification.tsv` | `results/qc/` | Classificação final, frase conclusiva e justificativa objetiva | É o artefato mais direto para decisão operacional de liberação |
| `tumor_purity_estimate.tsv` | `results/tables/` | Estimativa heurística de pureza tumoral, confiança, pico dominante de VAF e limitações | Use como contexto exploratório da amostra, nunca como medida isolada validada |
| `executive_summary.md` | `results/qc/` | Resumo textual curto da amostra | Útil para registro, auditoria e anexação em documentação interna |
| `panel_variant_status.tsv` | `results/tables/` | Relação entre genes do painel, cobertura e presença de variantes | Ajuda a responder quais genes estavam adequados e se tiveram ou não variantes observadas |
| `qc_report.html` | `results/reports/` | Relatório principal com resumo executivo, tabelas, alertas e gráficos | É o documento mais indicado para revisão humana e discussão técnico-laboratorial |
| `qc_report.md` | `results/reports/` | Versão em Markdown do relatório | Útil para versionamento, diff de texto e integração com sistemas simples |
| `plot_manifest.tsv` | `results/reports/` | Índice interno dos gráficos gerados no relatório | Serve mais para rastreabilidade e integração futura do que para revisão manual |

#### Arquivos de coorte

| Arquivo | Onde fica | O que contém | Como interpretar / quando usar |
| --- | --- | --- | --- |
| `cohort_sample_summary.tsv` / `cohort_sample_summary.csv` | `results/tables/` | Resumo consolidado de todas as amostras processadas no run | Melhor ponto de partida para comparar desempenho entre amostras e lotes |
| `cohort_gene_coverage_heatmap.html` | `results/reports/` | Heatmap gene x amostra com cobertura média | Útil para detectar genes sistematicamente fracos, problemas por lote e padrões recorrentes do painel |
| `tool_versions.tsv` | `results/logs/` ou `results/qc/` | Registro de versões das ferramentas usadas no run | Importante para rastreabilidade, validação e reprodutibilidade laboratorial |
| `trace.txt`, `execution_report.html`, `timeline.html`, `pipeline_dag.html` | `results/logs/` | Metadados operacionais de execução do Nextflow | Use para troubleshooting, auditoria computacional e análise de performance |

## Exemplo de output esperado

```text
results/
├── coverage/
│   ├── sample1.coverage_summary.tsv
│   ├── sample1.gene_coverage.tsv
│   ├── sample1.interval_coverage.tsv
│   ├── sample1.low_coverage_genes.tsv
│   └── sample1.low_coverage_intervals.tsv
├── qc/
│   ├── qc_rules.json
│   ├── sample1.classification.tsv
│   ├── sample1.executive_summary.md
│   └── sample1.qc_flags.tsv
├── reports/
│   ├── cohort_gene_coverage_heatmap.html
│   ├── sample1.plot_manifest.tsv
│   ├── sample1.qc_report.html
│   └── sample1.qc_report.md
├── tables/
│   ├── cohort_sample_summary.tsv
│   ├── panel_genes.tsv
│   ├── sample1.panel_variant_status.tsv
│   ├── sample1.sample_summary.tsv
│   └── sample1.tumor_purity_estimate.tsv
└── variants/
    ├── sample1.variant_gene_counts.tsv
    ├── sample1.vcf_summary.tsv
    ├── sample1.variants.tsv
    ├── sample1.variants_nonpass.tsv
    └── sample1.variants_pass.tsv
```

## DAG simplificado

Ver `docs/dag.txt`.

## Estimativa heurística de pureza tumoral

O pipeline inclui um módulo opcional de estimativa exploratória de pureza tumoral usando apenas os arquivos já disponíveis no fluxo, isto é, `BAM + VCF + BED` e as métricas derivadas desses artefatos. No estado atual do projeto, ele vem habilitado por padrão e gera uma tabela própria por amostra.

### Objetivo

Essa estimativa não substitui métodos dedicados de pureza/ploidia, mas pode fornecer um valor aproximado para contextualização técnica da amostra, especialmente em painéis somáticos sem normal pareado.

### Estratégia proposta

A heurística parte do princípio de que variantes somáticas clonais heterozigóticas em regiões diploides tendem a apresentar `VAF` proporcional à pureza tumoral. Sob essa hipótese simplificada:

- `pureza tumoral aproximada ~ 2 x VAF clonal dominante`
- exemplo: um pico clonal em `VAF ~ 0.22` sugere pureza tumoral aproximada de `44%`

### Seleção das variantes usadas na estimativa

A estimativa usa preferencialmente variantes que atendam aos critérios abaixo:

- `FILTER=PASS`
- `DP` alto, com threshold configurável
- `VAF` em faixa plausível para variante somática clonal
- variante localizada dentro das regiões alvo do painel
- preferência por `SNVs`, evitando misturar indiscriminadamente com `INDELs`
- exclusão de variantes com suporte técnico fraco, strand bias alto ou baixa confiabilidade

Filtros práticos desta primeira versão:

- `DP >= 200`
- `0.05 <= VAF <= 0.45`
- usar apenas `SNVs PASS`
- excluir variantes marcadas como `weak_support=true`
- excluir variantes `outside_panel=true`

Esses limites devem ser tratados como parâmetros configuráveis e ajustados após validação local do laboratório.

### Como a pureza é estimada

Fluxo conceitual:

1. selecionar variantes elegíveis
2. construir a distribuição de `VAF` dessas variantes
3. identificar um pico clonal dominante
4. calcular `pureza_aproximada = min(1.0, 2 x VAF_pico)`
5. classificar a confiança da estimativa

Sugestão inicial de classes de confiança:

- `NOT_ESTIMABLE`: poucas variantes elegíveis ou distribuição muito dispersa
- `LOW_CONFIDENCE`: sinal fraco, poucos eventos ou pico mal definido
- `MODERATE_CONFIDENCE`: número razoável de variantes elegíveis e pico relativamente estável

### Saídas atuais do módulo

A implementação atual gera:

- `tables/<sample>.tumor_purity_estimate.tsv`
- incorporação do resultado resumido em `tables/<sample>.sample_summary.tsv` e `tables/<sample>.sample_summary.csv`
- incorporação da estimativa no `reports/<sample>.qc_report.html` e `reports/<sample>.qc_report.md`

Campos atuais na tabela:

- `sample_id`
- `eligible_variants`
- `dominant_vaf_peak`
- `estimated_tumor_purity`
- `confidence_class`
- `method`
- `assumptions`
- `limitations`

### Parâmetros relacionados

- `enable_tumor_purity`
- `tumor_purity_min_dp`
- `tumor_purity_min_vaf`
- `tumor_purity_max_vaf`
- `tumor_purity_min_variants`
- `tumor_purity_bin_width`

### Interpretação recomendada

A interpretação laboratorial deve deixar claro que esse valor serve para contexto técnico e não deve ser usado isoladamente para decisão clínica ou molecular sem validação adicional.

### Limitações específicas dessa abordagem

As limitações abaixo devem constar idealmente no relatório ou outros documentos disponibilizados posteriormente quando o módulo for ativado:

- a abordagem assume, de forma simplificada, variantes clonais heterozigóticas em regiões diploides
- alterações de número de cópias, LOH e aneuploidia podem distorcer fortemente o `VAF` e enviesar a pureza estimada
- sem amostra normal pareada, a separação entre variante somática e germinativa é imperfeita
- painéis alvo têm cobertura genômica restrita e podem não conter variantes suficientes para uma inferência robusta
- subclonalidade tumoral pode gerar múltiplos picos de `VAF`, dificultando a identificação do componente clonal dominante
- `INDELs`, hotspots específicos e variantes em regiões tecnicamente difíceis podem produzir distribuições artificiais de `VAF`
- filtros excessivamente permissivos aumentam ruído; filtros excessivamente rígidos podem deixar a amostra sem variantes elegíveis
- o valor obtido deve ser tratado como estimativa aproximada, não como medida validada de pureza tumoral
- a estimativa não substitui abordagens dedicadas como modelos com `CNV`, `BAF`, `tumor-normal`, `CNVkit` ou `PureCN`

### Recomendação de uso

A recomendação pragmática é introduzir esse módulo como funcionalidade opcional, explicitamente marcada como:

- `estimativa heurística`
- `uso exploratório`
- `não validada para decisão isolada`

Em contexto laboratorial real, o ideal é validar essa estimativa contra casos com pureza conhecida por patologia, revisão morfológica, métodos ortogonais ou ferramentas dedicadas de pureza/ploidia.

## Limitações documentadas

- Se o VCF tiver múltiplas amostras, o parser usa a primeira amostra do header
- Campos `AF`, `DP`, `MQ` e `SB` são extraídos quando disponíveis; se ausentes, o pipeline registra vazio e segue
- Sobreposição complexa de intervalos no BED é suportada por busca local, adequada para painéis direcionados usuais
- O modo mínimo requer `BAM + BAI + VCF + BED`; FASTA é opcional e só é necessária para gerar dados sintéticos do `test_data`

## Implementação para próxima versão

- múltiplas classes de genes críticos com thresholds específicos
- blacklist de regiões sistematicamente difíceis
- análise específica de hotspots
- regras UMI-aware
- integração com anotação terciária
- suporte explícito a tumor-normal
- incorporação opcional de MultiQC
- priorização por regiões clinicamente mandatórias

## Comando mínimo de teste

```bash
cd test_data
bash generate_inputs.sh
cd ..
nextflow run main.nf -profile standard --input_dir test_data/input --outdir results_test
```
