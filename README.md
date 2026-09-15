# IBD Metagenome Strain-Level Diversity Analysis

This repository contains R scripts used to analyze strain-level genetic
variation across 20 microbial species, comparing gut metagenomes from IBD
patients (Crohn's Disease [CD] and Ulcerative Colitis [UC]) to healthy
Control metagenomes. Input data are per-species outputs from
[inStrain](https://instrain.readthedocs.io/), summarizing gene- and
variant-level metrics from metagenomic read alignments.

Each script was originally run **once per species**, as a job in a SLURM array — the
array wrapper sets a `file` variable to the path of that species' input
table before sourcing/running the script.

## Two central metrics

- **ND (nucleotide diversity)** — the `nucl_diversity` column from
  inStrain's gene-level output. A measure of within-sample genetic
  diversity at a gene or genome.
- **NS (nonsynonymous-variant enrichment)** — derived from the counts of
  nonsynonymous vs. synonymous variants (`SNV_N_count`/`SNV_S_count` for
  segregating variants, `SNS_N_count`/`SNS_S_count` for fixed
  substitutions), and the corresponding ratios `pNpS_variants` and
  `dNdS_substitutions`. Reports whether variation in a gene skews toward
  protein-changing (nonsynonymous) mutations, a signal of selection
  pressure.

## Shared input: the gene_info table

Scripts `2`,`9` (ND & NS), and `10` all read an inStrain gene-level
table with no header, so the column names are declared manually in each
script:

```
filename_scaffold, gene, gene_length, coverage, breadth, breadth_minCov,
nucl_diversity, start, end, direction, partial, dNdS_substitutions,
pNpS_variants, SNV_count, SNV_S_count, SNV_N_count, SNS_count,
SNS_S_count, SNS_N_count, divergent_site_count
```

## Objects expected in the environment

None of these scripts are self-contained — they assume certain objects
have already been created (typically by a wrapper script that sources them
inside the SLURM array job, after setting up sample metadata and paths):

- `file` — path to the current species' input table (set per array task)
- `genome_presence` - data frame with the qc'ed species-sample pairs 
- `samples_phenotype` — data frame mapping sample `label` to `subject`,
  `Phenotype` (`Control`/`CD`/`UC`), and `cohort` (`Pediatric`/`Adult`)
- `output_wd` — base directory for writing results
- `model_name` (scripts `2`, `10`) or `task` (scripts `4`, `9` ND/NS) —
  subfolder under `output_wd` for that script's output files

## Common per-species preprocessing (repeated across scripts)

Each script that reads the gene_info table performs a similar cleanup and
QC pipeline before its specific analysis:
1. Parse the sample label and genome ID out of `filename_scaffold`
   (with some fixes for a specific genome ID format in the UHGG v1
   catalogue).
2. Filter to genes covered at breadth ≥ 0.5.
3. Fill in occasional missing `gene_length` values.
4. Join in sample Phenotype/cohort metadata.
5. Where multiple samples exist per subject (or multiple gene instances
   per subject), keep only the most complete one (highest `breadth_minCov`).
6. Restrict to genome-sample pairs confirmed present via the
   `genome_presence` table.
7. Apply a `breadth_minCov` threshold (typically 0.2) to exclude poorly
   covered gene/genome instances.

## What do these scripts do?

Scripts are numbered according to their place in the overall analysis
pipeline; not every numbered step is included here (e.g. earlier steps
that build the `samples_phenotype`, `genome_presence`, and upstream gene
z-score/candidate-gene lists are not part of this set).

| Script | Analysis stage | What it does |
|---|---|---|
| `2_run_gene_Weighted_avg_genome_sample_variance_model.R` | Genome-level ND model | Collapses per-gene ND into one weighted genome-level ND value per sample, then fits a linear model of ND ~ Phenotype (Control/CD/UC), adjusting for coverage, breadth, and cohort. Outputs model coefficients, adjusted group means, and a variance-partitioning breakdown. |
| `9_run_ND_genes_stats.R` | Gene-level ND follow-up | Extracts per-sample, per-gene diversity/enrichment metrics for genes previously flagged as ND-decreased ("purged") or ND-increased ("diversified") in IBD vs. Control. |
| `9_run_NS_genes_stats.R` | Gene-level NS follow-up | Extracts per-sample, per-gene metrics for genes previously flagged as NS-increased or NS-decreased in IBD vs. Control, and correlates nonsynonymous enrichment (pN/pS) with clinical inflammation severity (fecal calprotectin) per Phenotype and gene category. |
| `10_run_gene_SNV_NS_enrich_model.R` | Gene-level NS enrichment model | Binomial mixed model (`glmmTMB`) testing whether nonsynonymous SNV enrichment differs by Phenotype, with per-gene random effects. Produces both an overall (population-level) effect and per-gene z-scores used to define the NS-increased/decreased gene lists consumed by `9_run_NS_genes_stats.R`. |

**Typical order of use:** `2` (genome-level ND) and `10` (gene-level NS
model) produce the statistical results (and, for `10`, the candidate gene
lists) that feed into the gene-level follow-up scripts `9_run_ND_genes_stats.R`,
`9_run_NS_genes_stats.R`, and `4_run_snv_stats.R`.

## Notes on script `9_run_NS_genes_stats.R`

Two of the three result tables this script computes (mean gene counts
retained before/after QC filtering, and count-weighted mean pN/pS & dN/dS
per gene and Phenotype) have their `write_csv()` calls commented out, so
by default only the severity-correlation table
(`<species>_pnps_corr_severity.csv`) is written to disk.
