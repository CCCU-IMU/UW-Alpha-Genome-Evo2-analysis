# Methods summary for the UW selected12 model-based prioritization workflow

This file records the methodological assumptions needed to understand and reuse
the public code package. It is not a replacement for the full manuscript
Methods.

## Candidate variants

Candidate variants were prioritized inside pre-defined Ujimqin White cattle
selection-signal loci centered on or near `ASIP`, `PMEL`, `SMG6`, `TRIM59`,
`PLAG1` and `GDF11`. The private project used a candidate31 set and then chose
two representative SNPs per target gene/locus for the selected12 analysis.

The public repository includes only small example tables. Users should provide
their own candidate table if rerunning the workflow.

## Reference sequence and REF/ALT construction

In the private analysis, model input sequences were generated from the cattle
ARS-UCD2.0 reference genome. REF and ALT sequences were constructed with a
1,048,576 bp window, placing the SNP at 1-based offset 524,288 for non-edge
variants. The ALT sequence differed from the REF sequence only at the candidate
SNP base.

The public code does not distribute the reference FASTA or generated model input
sequences. Users should provide their own reference FASTA and generated sequence
summaries.

## Evo2 sequence-level evidence

Evo2 was used to derive sequence-level variant-effect evidence as a delta
likelihood:

```text
delta_likelihood = ALT_score - REF_score
abs_delta_likelihood = abs(delta_likelihood)
```

The private selected12 workflow used a context length of 16,384 bp for the
candidate31/selected12 summarization. The public workflow expects summarized
Evo2 fields such as `evo2_abs_delta_likelihood`; it does not distribute Evo2
weights or GPU job scripts.

## AlphaGenome regulatory-track evidence

AlphaGenome DNASE, ATAC and expression-related outputs were used as model-based
regulatory predictions. Expression-related outputs included CAGE, PROCAP and
RNA-seq prediction summaries.

Important interpretation note: AlphaGenome outputs are model predictions under a
human regulatory-track framework and are not measured cattle RNA or chromatin
data.

Variant effects were summarized as:

```text
delta = ALT_prediction - REF_prediction
center_abs_delta_sum = sum(abs(delta)) in SNP-centered window
center_max_abs_delta = max(abs(delta)) in SNP-centered window
```

The private analysis used a SNP-centered window of ±128 bp for the center
metrics. Large raw prediction arrays are not part of this public repository.

## UW allele-direction analysis

For each candidate SNP, REF/ALT allele frequencies were compared between UW and
control populations. A simplified public implementation uses group-level AF
fields and a configurable margin.

Direction-aware interpretation follows:

- UW-enriched ALT: use model `ALT - REF` directly.
- UW-enriched REF: use `REF - ALT`, i.e. multiply `ALT - REF` by `-1`.
- intermediate / low-confidence: do not assign a strong direction-aware effect.

## LD/haplotype follow-up

Some representative SNPs may not be strongly UW-enriched by themselves. The
private follow-up checked whether such SNPs were linked to regional top
delta-AF SNPs or shared the same regulatory window. The public script computes a
simple dosage-based LD r² for representative/top-variant pairs.

The working LD support threshold used in the manuscript workflow was:

```text
r² >= 0.30
```

## Targeted motif gain/loss

Targeted motif analysis compares local REF and ALT sequences around candidate
SNPs. The private workflow used curated transcription-factor motifs relevant to
pigmentation, development and stress/adaptation biology. Motif score changes
were interpreted as hypotheses only.

The public example script implements a minimal PWM scanner for transparent
demonstration. Users should document the motif database and version if applying
the workflow to real data, for example a specific JASPAR release.

## Normalization and integrated score

The public workflow follows the figure-level normalization used in the private
analysis:

- selection priority: min-max normalization
- Evo2 absolute delta likelihood: min-max normalization
- DNASE, ATAC, expression-related and motif metrics: `log1p(raw_abs_value)`
  followed by min-max normalization

Default weights:

```text
Evo2                 0.40
DNASE                0.20
ATAC                 0.16
Expression-related   0.06
Selection priority   0.10
Motif proxy          0.08
```

The integrated score is a ranking and visualization support score. It should not
be described as a raw biological effect size.

## Recommended manuscript wording

Use wording such as:

- `model-based prioritization`
- `candidate mechanism hypothesis`
- `predicted regulatory effect`
- `sequence-level likelihood effect`

Avoid wording that implies experimental validation, direct TF binding, measured
UW expression, or causal proof in cattle.
