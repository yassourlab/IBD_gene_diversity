##############################################################################
# 10_run_gene_SNV_NS_enrich_model.R
#
# PURPOSE
# -------
# For a single species, test whether nonsynonymous SNV enrichment (NS = the
# balance of nonsynonymous vs synonymous SNVs, i.e. SNV_N_count vs
# SNV_S_count) differs between IBD Phenotypes (CD, UC) and Control, across
# a curated set of genes, using a binomial mixed-effects model. Each gene
# gets its own random intercept and random Phenotype slope, so the model
# both estimates an overall (population-level) NS-enrichment effect of
# Phenotype and, from the per-gene random effects, identifies individual
# genes with the strongest Phenotype-associated shift toward more/fewer
# nonsynonymous variants (used to build the "NS_increased_genes" /
# "NS_decreased_genes" gene lists consumed by 9_run_NS_genes_stats.R).
#
# This script is run once per species as a SLURM array job: "file" is the
# path to that species' per-gene inStrain output (.tsv).
#
# EXPECTED OBJECTS IN ENVIRONMENT (set by the calling/wrapper script, not
# defined here):
#   file              - path to the species' gene-level inStrain .tsv
#   samples_phenotype - sample metadata (label -> subject, Phenotype, cohort)
#   output_wd         - base output directory
#   model_name        - subfolder name under output_wd for this model's
#                       results
#
# MAIN OUTPUTS (written to output_wd/model_name/):
#   <species>_model_res.csv            - fixed-effect coefficients from the
#                                         binomial GLMM (overall Phenotype
#                                         effect on nonsyn/syn SNV ratio)
#   <species>_random_gene_estimates.csv - per-gene random-effect estimates,
#                                         standard errors, and z-scores for
#                                         the CD and UC Phenotype effects
#                                         (used to flag individual
#                                         NS-enriched/depleted genes)
##############################################################################

# Minimum breadth at min-coverage required for a gene to be considered
# reliably genotyped in a sample.
breadth_minCov_TH <- 0.2

# Extract the species name from the input file path.
sp_name <- str_remove(str_split_fixed(file,"/",n=12)[,12],"\\.tsv")

# Genome presence table for this species: which genome was called present
# in which sample.
path="/sci/backup/morani/lab/Projects/IBD_diversity/UHGGv1_genomes/2.merged_instrain_output/genome/dfs_plots/genome_presence/"
genome_presence <- read.delim(file=paste0(path,sp_name,"_genome_sample_presence.csv"), header=T, sep=",") %>% dplyr::select(label,genome)

# Column names for inStrain's gene_info.tsv output (see script 2 for
# details on each column). Key columns for this script:
#   SNV_N_count / SNV_S_count -> nonsynonymous / synonymous SNV counts per
#                                 gene-sample, the response variable for
#                                 the NS-enrichment model below
column_names <- c("filename_scaffold", "gene","gene_length",
"coverage","breadth","breadth_minCov","nucl_diversity","start",
"end","direction","partial","dNdS_substitutions","pNpS_variants",
"SNV_count","SNV_S_count","SNV_N_count","SNS_count","SNS_S_count","SNS_N_count",
"divergent_site_count")


# Read and parse the raw per-gene inStrain output for this species (same
# parsing/genome-ID cleanup logic as in the other scripts).
sp_gene_info.raw <- read.delim(file=file, header = F, col.names = column_names) %>%
   mutate(gene_type = "non-marker",
          label = str_extract(filename_scaffold, "(?<=/)[^/]+(?=/output)"),
          genome = str_split_fixed(filename_scaffold, ":", n=2)[,2],
          genome = ifelse(str_detect(genome,"CAX"),"_964248265_",genome),
          genome = str_split_fixed(genome,"_", n=3)[,2], 
          genome = paste0("GUT_", genome),
          gene = str_replace(gene,"CAXTYD010000","GUT_964248265_")) %>% 
  dplyr::select(-filename_scaffold) %>% 
  filter(breadth >= 0.5)

#print(head(sp_gene_info.raw))
#print(nrow(sp_gene_info.raw))

# Build the per-subject, per-gene QC-passing gene table: fill in missing
# gene lengths, attach Phenotype/cohort metadata, keep the most complete
# instance of each gene per subject, restrict to genome-sample pairs with
# confirmed genome presence, and apply the breadth_minCov threshold.
sp_gene_info <- sp_gene_info.raw %>%
  # fill in the missing gene lengths (because later we use them as weights)
  group_by(gene) %>%
  fill(gene_length, .direction = "downup") %>% ungroup() %>%
  inner_join(samples_phenotype, by = "label") %>%
  # select the most complete sample per subject # breadth_min_cov #
  group_by(subject,gene) %>%
  arrange(desc(breadth_minCov)) %>%
  slice_head(n=1) %>% ungroup() %>% 
  mutate(Phenotype = factor(Phenotype , levels = c("Control", "CD", "UC")),
  cohort = factor(cohort, levels = c("Pediatric","Adult"))) %>%
  # select samples that have passed the 0.5 breadth threshold
  inner_join(
    dplyr::select(genome_presence, c("label","genome")),
  by = c("label","genome")) %>%
  filter(breadth_minCov >= breadth_minCov_TH )

# Restrict the model to a pre-selected list of genes: those retained by an
# earlier Bayesian gene-weighting step of the ND model (2-stage Bayesian
# model fit on Control samples), which serves here as a general-purpose
# "reliable gene" filter (genes with stable/well-supported estimates)
# rather than something specific to NS.
path="/sci/backup/morani/lab/Projects/IBD_diversity/UHGGv1_genomes/2.merged_instrain_output/gene/models/nucl_diversity/gene_weights_2stage_bayesian/"
rel_genes.list <- read.delim(file=paste0(path,"bayesian_model_Control_all.csv"), header=T, sep=",") %>% pull(gene) %>% unique()

# Build the model input: restrict to the reliable-gene list, impute missing
# dN/dS and pN/pS as 0 (genes with no substitutions/variants observed),
# collapse Phenotype into a binary Control vs IBD "status" (in addition to
# keeping the 3-level Phenotype), and z-scale log coverage as a covariate.
# Then, per gene, require at least 10 samples with any nonsyn/syn SNVs
# observed and at least 5 total nonsyn+syn SNVs across samples, so the
# per-gene random effects below are estimated from genes with enough
# variant evidence.
model_input <- sp_gene_info %>%
filter(gene %in% rel_genes.list) %>%
mutate(dNdS_substitutions = ifelse(is.na(dNdS_substitutions),0,dNdS_substitutions),
       pNpS_variants = ifelse(is.na(pNpS_variants),0,pNpS_variants),
       status = ifelse(Phenotype == "Control","Control","IBD"),
       status = factor(status , levels = c("Control", "IBD")),
       Phenotype = factor(Phenotype, levels = c("Control","CD","UC")),
       log_cov_scaled = as.vector(scale(log(coverage)) )) %>% 
# filtering out genes with little information
group_by(gene) %>%
  # n samples with more than 0
  filter(sum((SNV_N_count + SNV_S_count) > 0) >= 10) %>%
  # total sum of variants more than 5
  filter(sum(SNV_N_count + SNV_S_count) >= 5) %>% 
unique()

# Binomial GLMM: for each gene-sample, model the nonsynonymous vs
# synonymous SNV counts (cbind(SNV_N_count, SNV_S_count)) as a function of
# Phenotype, adjusting for dN/dS (fixed substitution signal), scaled log
# coverage, and breadth_minCov. Each gene gets its own random intercept and
# random Phenotype slope (1 + Phenotype || gene), letting the
# Phenotype-NS-enrichment relationship vary by gene. The "||" (rather than
# "|") drops estimation of the correlation between the random intercept and
# slopes, which was needed to avoid a non-positive-definite Hessian
# (model convergence issue) with the full covariance structure.
raw_counts.model_across_genes <-
    glmmTMB::glmmTMB(
      # double pipe || removes covariance estimation -> otherwise problem with non-positive-definite Hessian matrix
      cbind(SNV_N_count, SNV_S_count) ~ Phenotype + dNdS_substitutions + log_cov_scaled + breadth_minCov + (1 + Phenotype || gene),
      data = model_input,
      family = binomial(),
      #family = betabinomial(),
      control = glmmTMBControl( optimizer = nlminb)
)


# diagnose the fit
# Run glmmTMB's built-in diagnostics to flag convergence issues, singular
# fits, or other problems with this species' model before trusting its
# output.
cat("\nIs there anything wrong with this fit according to diagnose()?\n")
print(diagnose(raw_counts.model_across_genes))
   
# Tidy fixed-effect coefficients (overall Phenotype effect on NS
# enrichment, plus covariate effects), for the population-level result.
raw_counts.model_across_genes.df <- tidy(raw_counts.model_across_genes) %>% mutate(genome = sp_name)

model_res_path <- paste0(output_wd, "/", model_name, "/", sp_name, "_model_res.csv")
write_csv(raw_counts.model_across_genes.df, model_res_path, quote = "none")


# Extract the per-gene random effects (deviations from the overall
# Phenotype effect for each gene), along with their conditional variances,
# to identify individual genes with unusually strong NS enrichment shifts.
raw_counts.model_across_genes.re <- ranef(raw_counts.model_across_genes, condVar = TRUE)$cond$gene

# covariance matrix of the effects
# Conditional variance-covariance array for the random effects, one matrix
# slice per gene; used below to get each gene's standard error for the CD
# and UC random slopes.
condVar_cv <- attr(raw_counts.model_across_genes.re, "condVar")

# Standard errors of each gene's random CD-effect (row/col 2) and
# UC-effect (row/col 3) on the log-odds of nonsynonymous vs synonymous SNVs.
cd_se <- sapply(1:dim(condVar_cv)[3], function(i) {
  sqrt(condVar_cv[2,2,i])
})

uc_se <- sapply(1:dim(condVar_cv)[3], function(i) {
  sqrt(condVar_cv[3,3,i])
})

# estimates
# Per-gene random-effect point estimates for the CD and UC Phenotype
# slopes (i.e. how much more/less nonsynonymous-enriched that gene is in
# CD/UC relative to Control, beyond the population-average effect).
cd_effect <- raw_counts.model_across_genes.re$PhenotypeCD
uc_effect <- raw_counts.model_across_genes.re$PhenotypeUC

# standard errors
raw_counts.model_across_genes.re$cd_se <- cd_se
raw_counts.model_across_genes.re$uc_se <- uc_se

# z scores
# Per-gene z-scores (estimate / SE) for the CD and UC effects -- these are
# the values used downstream (e.g. in the wrapper script that builds
# NS_increased_genes.txt / NS_decreased_genes.txt) to flag genes with a
# significant, gene-specific shift in nonsynonymous enrichment.
raw_counts.model_across_genes.re$cd_z <- cd_effect / cd_se
raw_counts.model_across_genes.re$uc_z <- uc_effect / uc_se

raw_counts.model_across_genes.re$genome <- sp_name
raw_counts.model_across_genes.re$gene <- rownames(raw_counts.model_across_genes.re)


gene_res_path <- paste0(output_wd, "/", model_name, "/", sp_name, "_random_gene_estimates.csv")
write_csv(raw_counts.model_across_genes.re, gene_res_path, quote = "none")
