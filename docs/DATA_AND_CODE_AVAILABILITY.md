# Draft Data and Code Availability wording

The code used to organize the UW selected12 model-based candidate variant
prioritization workflow is available at the accompanying GitHub repository and
archived on Zenodo [repository URL and DOI to be added after release]. The
public repository contains lightweight scripts, documentation and example input
tables for integrating selection-priority, Evo2, AlphaGenome prediction-summary,
allele-direction, LD/haplotype and motif gain/loss evidence. Raw sequencing
data, private genotype files, full VCFs, AlphaGenome API keys, large raw
prediction arrays and cluster execution logs are not distributed in the code
repository. Users should provide their own reference genome FASTA, VCF-derived
allele-frequency summaries, Evo2 output summaries, AlphaGenome output summaries
and motif files to rerun the workflow. AlphaGenome outputs in this study are
model predictions under a human regulatory-track framework and are not measured
cattle RNA expression or chromatin-accessibility data.
