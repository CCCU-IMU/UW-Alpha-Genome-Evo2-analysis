# UW candidate variant prioritization workflow

This repository contains a minimal, public, reviewer-readable code package for
the **Ujimqin White cattle (UW) selected12 model-based candidate variant
prioritization workflow**.

The workflow integrates selection-priority evidence, curated candidate variants,
Evo2 sequence-level delta likelihood, AlphaGenome DNASE/ATAC/CAGE/PROCAP/RNA-seq
prediction summaries, UW allele-direction analysis, LD/haplotype follow-up,
targeted motif gain/loss and evidence-matrix plotting.

These analyses should be interpreted as **model-based prioritization and
candidate mechanism hypotheses**, not functional validation. AlphaGenome outputs
are model predictions under a human regulatory-track framework and are **not**
measured cattle RNA expression, DNASE, ATAC, CAGE or PROCAP data.

## What is included

- Small example TSV/FASTA files in `examples/`.
- Lightweight scripts in `scripts/` that reproduce the public table-level
  workflow.
- Shared utilities in `src/uw_variant_prioritization/`.
- A Methods-oriented summary in `METHODS_SUMMARY.md`.
- A code/data availability draft in `docs/DATA_AND_CODE_AVAILABILITY.md`.

## What is not included

This public repository intentionally does **not** include:

- Raw FASTQ, BAM, VCF or genotype matrices.
- Private UW sample metadata.
- AlphaGenome API keys or server logs.
- Large AlphaGenome raw NPZ arrays.
- Evo2 checkpoints or GPU queue submission details.
- Absolute server paths or private database paths.

Users should provide their own reference FASTA, VCF-derived allele-frequency
tables, Evo2 outputs, AlphaGenome output summaries and motif files when applying
the workflow to another dataset.

## Repository layout

```text
examples/
  selected12_variants.example.tsv
  candidate31_curated_variants.example.tsv
  sample_group.example.tsv
  integrated_evidence_table.example.tsv
  allele_frequency.example.tsv
  ld_pairs.example.tsv
  dosage.example.tsv
  ref_alt_sequences.example.fa
  motif_pwm.example.tsv
scripts/
  01_build_integrated_evidence.py
  02_allele_direction_from_af.py
  03_ld_followup_from_dosage.py
  04_targeted_motif_gain_loss.py
  05_plot_evidence_figures.py
src/uw_variant_prioritization/
  utils.py
METHODS_SUMMARY.md
environment.yml
```

## Quick start

Create the environment:

```bash
conda env create -f environment.yml
conda activate uw_variant_prioritization
export PYTHONPATH="$PWD/src:${PYTHONPATH}"
```

Run the example workflow:

```bash
python scripts/01_build_integrated_evidence.py \
  --input examples/integrated_evidence_table.example.tsv \
  --output results/integrated_evidence.ranked.tsv

python scripts/02_allele_direction_from_af.py \
  --input examples/allele_frequency.example.tsv \
  --output results/allele_direction.tsv

python scripts/03_ld_followup_from_dosage.py \
  --pairs examples/ld_pairs.example.tsv \
  --dosage examples/dosage.example.tsv \
  --output results/representative_vs_top_deltaAF_LD.tsv

python scripts/04_targeted_motif_gain_loss.py \
  --fasta examples/ref_alt_sequences.example.fa \
  --motifs examples/motif_pwm.example.tsv \
  --output results/targeted_motif_gain_loss.tsv

python scripts/05_plot_evidence_figures.py \
  --input examples/integrated_evidence_table.example.tsv \
  --output-prefix figures/example_evidence
```

Each script also supports `--help` and `--dry-run`.

## Workflow overview

1. **Candidate definition.** Candidate variants are defined inside prior
   selection-signal loci and summarized as a candidate31 table and a selected12
   representative set.
2. **Sequence-level model evidence.** Evo2 delta likelihood is used as
   sequence-level support. The public script expects an already summarized
   `evo2_abs_delta_likelihood` field.
3. **Regulatory-track model evidence.** AlphaGenome DNASE, ATAC and
   expression-related outputs are summarized as center-window ALT-minus-REF
   scores. The public workflow uses table summaries only.
4. **Allele direction.** UW allele-frequency summaries are used to determine
   whether the UW-biased allele is REF, ALT, intermediate or low-confidence.
5. **LD/haplotype follow-up.** Representative SNPs can be linked to regional
   top delta-AF variants using dosage-based LD r².
6. **Motif gain/loss.** REF and ALT local sequences are scanned against
   targeted motifs to generate motif gain/loss hypotheses.
7. **Evidence integration.** Long-tailed regulatory metrics are log1p
   transformed, min-max normalized and combined with fixed weights for a
   ranking-oriented integrated score.

## Required input columns

The integrated evidence script expects a TSV with at least:

- `variant_id`, `candidate_gene`, `chrom`, `pos`, `ref`, `alt`
- `priority_score`
- `evo2_abs_delta_likelihood`
- `dnase_center_abs_delta_sum`
- `atac_center_abs_delta_sum`
- `expression_related_abs_score`
- `targeted_motif_abs_normalized_delta`

The allele-direction script expects:

- `variant_id`, `ref`, `alt`
- `UW_AF`
- one or more control population AF fields such as `Angus_AF`,
  `Charolais_AF`, `MoOD_AF`, `MoSN_AF`, `Mongolia_AF` and `control_AF`.

## Interpretation guardrails

- Evo2 is treated as sequence-level likelihood evidence.
- AlphaGenome outputs are treated as regulatory-track predictions under a human
  regulatory-track framework applied to cattle sequence inputs.
- Motif gain/loss is an exploratory transcription-factor hypothesis.
- UW allele direction converts ALT-REF model output into a UW-biased allele
  direction only when allele-direction support is clear.
- None of these computational analyses is functional validation in cattle.

## License

Code is released under the MIT License. Example data are provided only to
demonstrate table formats and should not be interpreted as a public release of
the underlying private genotype dataset.
