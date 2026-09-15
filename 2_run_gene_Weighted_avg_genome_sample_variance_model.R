##############################################################################
# 2_run_gene_Weighted_avg_genome_sample_variance_model.R
#
# PURPOSE
# -------
# For a single bacterial species, collapse per-gene nucleotide diversity (ND)
# values (inStrain gene_info output) into one ND value per genome-sample,
# using a coverage- and length-weighted average across genes. This weighted
# genome-level ND is then modeled against IBD Phenotype (Control/CD/UC),
# while adjusting for sequencing coverage, breadth, and cohort (Pediatric/
# Adult), to test whether within-species nucleotide diversity differs
# between healthy and IBD metagenomes.
#
# This script is run once per species as a SLURM array job: "file" is the
# path to that species' per-gene inStrain output (.tsv), passed in from the
# array wrapper.
#
# EXPECTED OBJECTS IN ENVIRONMENT (set by the calling/wrapper script, not
# defined here):
#   file              - path to the species' gene-level inStrain .tsv (SLURM
#                        array input)
#   genome_presence   - data frame containing qc'ed species-sample pairs
#   gp_path           - path to the genome_presence file
#   samples_phenotype - data frame mapping sample "label" -> "subject",
#                        "Phenotype" (Control/CD/UC), "cohort"
#                        (Pediatric/Adult)
#   output_wd         - base output directory
#   model_name         - subfolder name under output_wd for this model's
#                        results
#
# MAIN OUTPUTS (written to output_wd/model_name/):
#   <species>_genome_nd_w_wo_weights.csv   - per-sample comparison of
#                                             weighted vs. unweighted ND/
#                                             coverage/breadth (QC check)
#   <species>_input_and_residuals.csv      - model input data with fitted
#                                             residuals and adjusted ND
#   <species>_adj_means.csv                - estimated marginal means of ND
#                                             per Phenotype (from emmeans)
#   <species>_weighted_gene_model_res.csv  - linear model coefficients,
#                                             sample sizes, and variance
#                                             explained per term
##############################################################################

# Minimum fraction of the genome that must be covered (at the minimum
# coverage threshold used by inStrain, breadth_minCov) for a genome-sample
# pair to be considered reliably present/detected.
calc_genome_breadth_minCov_TH <- 0.2


# Extract the species name from the input file path
# (file naming convention: .../<species_name>.tsv)
sp_name <- str_remove(str_split_fixed(file,"/",n=12)[,12],"\\.tsv")

# genome presence
# Table of which genome (strain/representative) was called present in which
# sample for this species, produced upstream in the pipeline. Used below to
# restrict the analysis to genome-sample pairs that actually passed the
# genome-level presence/detection filter.
genome_presence <- read.delim(file=paste0(gp_path,sp_name,"_genome_sample_presence.csv"), header=T, sep=",") %>% dplyr::select(label,genome)

print(head(samples_phenotype))

# gene ND from the SPA run
#markers_info <- read.delim(file = "/sci/backup/morani/lab/Projects/IBD_diversity/markers_gene_info.csv", header =T, sep=",") %>% rename(marker_id = gene)
#print(head(markers_info))

# Column names for inStrain's gene_info.tsv output (this pipeline's tables
# have no header row, so the names are supplied manually here and reused
# across scripts).
# Key columns for this script:
#   nucl_diversity   -> per-gene nucleotide diversity (ND)
#   coverage         -> per-gene mean read coverage
#   breadth_minCov   -> fraction of the gene covered at the min-coverage
#                       threshold (used as a per-gene reliability weight)
#   gene_length      -> gene length in bp (used as a weight, alongside
#                       breadth_minCov, so that longer/more-covered genes
#                       contribute more to the genome-level average)
column_names <- c("filename_scaffold", "gene","gene_length",
"coverage","breadth","breadth_minCov","nucl_diversity","start",
"end","direction","partial","dNdS_substitutions","pNpS_variants",
"SNV_count","SNV_S_count","SNV_N_count","SNS_count","SNS_S_count","SNS_N_count",
"divergent_site_count")


# Read and parse the raw per-gene inStrain output for this species.
sp_gene_info.raw <- read.delim(file=file, header = F, col.names = column_names) %>%
   mutate(gene_type = "non-marker",
          # sample label is embedded in the scaffold/file path
          label = str_extract(filename_scaffold, "(?<=/)[^/]+(?=/output)"),
          # genome (representative genome ID) is embedded after the ":" in
          # the scaffold name; the block below normalizes a few genome ID
          # formatting quirks specific to this genome catalogue (UHGG v1)
          genome = str_split_fixed(filename_scaffold, ":", n=2)[,2],
          genome = ifelse(str_detect(genome,"CAX"),"_964248265_",genome),
          genome = str_split_fixed(genome,"_", n=3)[,2], 
          genome = paste0("GUT_", genome),
          gene = str_replace(gene,"CAXTYD010000","GUT_964248265_")) %>%
  dplyr::select(-filename_scaffold) %>%
  # remove genes that are not covered by at least half
  filter(breadth >= 0.5) %>%
  # fill in the missing gene lengths (do not know the reason)
  group_by(gene) %>%
  fill(gene_length, .direction = "downup") %>% ungroup() %>%
  dplyr::select(gene, gene_length, genome, label, nucl_diversity, coverage, breadth_minCov) %>%
  # remove NAs because they will disrupt the mean and the weighing
  na.omit()

print(head(sp_gene_info.raw))

# Per-sample "total weight" (sum across all genes of breadth_minCov *
# gene_length), i.e. roughly the total number of reliably-covered gene
# sites in that sample. This is later used as the observation-level weight
# in the linear model (samples with more covered sequence get more weight).
sample_weights <- sp_gene_info.raw %>%
group_by(label) %>%
reframe(weight = sum(breadth_minCov * gene_length))

# select only relevant gene/genome-sample pairs 
# Keep only gene-level records for (sample, genome) pairs where that genome
# was actually called present in that sample (per genome_presence table),
# and attach each sample's Phenotype/cohort metadata.
rel_gene_genomes <- sp_gene_info.raw %>%
inner_join(samples_phenotype, by = "label") %>%
# select samples that have passed the 0.5 breadth threshold
inner_join(
    dplyr::select(genome_presence, c("label","genome")),
   by = c("label","genome"))

print(head(rel_gene_genomes))

#stop()

# --- UNWEIGHTED genome-level metric (for comparison only) ---
# Simple, unweighted per-sample mean of nucleatide diversity/coverage/
# breadth_minCov across the genes retained above. This is *not* used in the
# model; it's computed purely to sanity-check the weighted approach below.
sp_gene_info.unW <- sp_gene_info.raw %>%
 inner_join(
   dplyr::select(rel_gene_genomes, c("label","gene")) ,
  by = c("label" ,"gene")) %>% 
  # recalculate nucl_diversity and coverage per sample using only marker genes
  group_by(label) %>%
  reframe(calc_genome_metric = mean(nucl_diversity),
          calc_genome_coverage = log(mean(coverage)),
          calc_genome_breadth_minCov = mean(breadth_minCov))
  

# --- WEIGHTED genome-level metric (the one actually used downstream) ---
# For each sample, collapse all genes into one genome-level nucleotide
# diversity value, weighting each gene's contribution by
# (breadth_minCov * gene_length), i.e. by how many reliably-covered sites
# that gene contributes. Coverage is also weighted and log-transformed
# (log coverage is used as a covariate to correct for sequencing depth
# effects on diversity estimates).
sp_gene_info.W <- sp_gene_info.raw %>%
  group_by(label) %>%
  # breadth_minCov * gene_length = number of gene sites used
  reframe(calc_genome_metric = sum(nucl_diversity * breadth_minCov * gene_length) / sum(breadth_minCov * gene_length),
          log_calc_genome_coverage = log( sum(coverage * breadth_minCov * gene_length) / sum(breadth_minCov * gene_length)),
          calc_genome_breadth_minCov = sum(breadth_minCov * gene_length) / sum(gene_length)) %>%
  mutate(genome = sp_name)

# show how weighing is better
# Merge the weighted and unweighted per-sample estimates side by side, so
# the effect of weighting can be inspected/plotted downstream (QC output).
wei_unwei.comparison <- sp_gene_info.unW %>%
inner_join(sp_gene_info.W, by = c("label"), suffix= c(".unweighted",".weighted"))

write_csv(wei_unwei.comparison, file = paste0(output_wd, model_name, "/", sp_name, "_genome_nd_w_wo_weights.csv"),
      quote = "none", num_threads = 2)


# --- Build the model input table ---
sp_gene_info <- sp_gene_info.W %>%
  inner_join(samples_phenotype, by = "label") %>%
  # remove very poorly represented genomes (< 0.2)
  filter(calc_genome_breadth_minCov >= calc_genome_breadth_minCov_TH) %>%
  # select the most complete sample per subject # breadth_min_cov #
  group_by(subject) %>%
  arrange(desc(calc_genome_breadth_minCov)) %>%
  slice_head(n=1) %>% ungroup() %>%
  inner_join(sample_weights, by = "label") %>%
  mutate(Phenotype = factor(Phenotype , levels = c("Control", "CD", "UC")),
  cohort = factor(cohort, levels = c("Pediatric","Adult"))) 

print(head(sp_gene_info))

# Model: weighted genome-level ND explained by IBD Phenotype, adjusting for
# log sequencing coverage, genome breadth (detection completeness), and
# cohort (Pediatric vs Adult), which can otherwise confound diversity
# comparisons.
formula <- as.formula( "calc_genome_metric ~ Phenotype + log_calc_genome_coverage + calc_genome_breadth_minCov + cohort" )

# input
model_input <- sp_gene_info %>% 
  dplyr::select(calc_genome_metric,label,Phenotype,log_calc_genome_coverage,calc_genome_breadth_minCov,cohort,weight) %>% 
  na.omit()

print(head(model_input))


# model
# Weighted least-squares linear model: each sample is weighted by the total
# amount of reliably-covered gene sequence it contributed (sample_weights),
# so samples with more/better-covered genes count more.
fit <- lm(formula, data = model_input, weights = model_input$weight)

### USE FIT
# compute mean of the covariates to predict adjusted nucl_diversity
# Build a version of the data where the nuisance covariates (coverage,
# breadth) are fixed at their sample means, so that predictions from the
# model isolate the Phenotype effect from coverage/breadth differences.
df_at_means <- model_input
df_at_means$log_calc_genome_coverage <- mean(model_input$log_calc_genome_coverage)
df_at_means$calc_genome_breadth_minCov <- mean(model_input$calc_genome_breadth_minCov)

# Calculate the adjusted value for every row -> old and new estimate for all samples, all fits
# Covariate-adjusted ND per sample = model residual + prediction at mean
# covariate values. This gives an ND value "as if" every sample had average
# coverage/breadth, useful for plotting/visual comparison across Phenotypes.
model_input$nucl_diversity_adjusted <- residuals(fit) + predict(fit, newdata = df_at_means)
model_input$genome <- sp_name

#Get the group means and confidence intervals -> 3 per fit
# Estimated marginal (least-squares) means of ND per Phenotype, adjusted for
# the other covariates in the model.
adj_means_Phenoype.df <- as.data.frame(emmeans(fit, ~ Phenotype))
adj_means_Phenoype.df$genome <- sp_name

# save
model_input_path <- paste0(output_wd, "/", model_name, "/", sp_name, "_input_and_residuals.csv")
write_csv(model_input, model_input_path, quote = "none")

# save
adj_means_Phenoype_path <- paste0(output_wd, "/", model_name, "/", sp_name, "_adj_means.csv")
write_csv(adj_means_Phenoype.df, adj_means_Phenoype_path, quote = "none")



# calculate var prop with other package
# Relative importance analysis (LMG method): partitions the model's R^2
# among the fixed-effect terms, to see how much of the explained variance in
# ND is attributable to Phenotype vs. the nuisance covariates.
var_res <- relaimpo::calc.relimp(fit, type = "lmg", rela = TRUE)

model.var <- data.frame(var_res$lmg) %>%
      rownames_to_column("term") %>% 
      rename(prop_var_fixed.comp = var_res.lmg) %>% 
      mutate(prop_var_fixed.comp = prop_var_fixed.comp*100, total_var = var_res$var.y,
      # missing columns in this package
      prop_var_fixed=NA,prop_var_random=NA, var.residual=NA, var.distribution=NA, var.dispersion=NA, var.intercept=NA) %>%
      mutate(term = str_replace(term,"cohort","cohortAdult"))

print(model.var)


# mean_n_samples
# Number of samples per Phenotype group actually used in the model (for
# reporting alongside effect sizes).
model_input.sample_stats <- model_input %>% group_by(Phenotype) %>% tally()

# combine output 
# Tidy up the fixed-effect coefficients (Phenotype effect estimates,
# standard errors, p-values), attach sample sizes, variance-explained
# figures, and the model's marginal R^2 (via MuMIn), into one results table
# for this species.
res_df <- tidy(fit, effects = "fixed") %>%
  mutate( genome = sp_name,
          term = str_replace(term, "Phenotype",""),
          term = str_replace(term, "\\(Intercept\\)","Control"),
          r2 = MuMIn::r.squaredGLMM(fit)[2]) %>% 
  left_join(model_input.sample_stats, by = c("term" = "Phenotype")) %>%
  full_join(model.var, by = "term") %>% 
  dplyr::select(genome,term,estimate,std.error,statistic,p.value,r2,n,
                total_var,prop_var_fixed.comp)

print(res_df)


# save
write_csv(res_df, file = paste0(output_wd, model_name, "/", sp_name, "_weighted_gene_model_res.csv"), quote = "none", num_threads = 2)

