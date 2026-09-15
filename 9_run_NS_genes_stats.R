##############################################################################
# 9_run_NS_genes_stats.R
#
# PURPOSE
# -------
# For a single species, pull out the per-sample, per-gene inStrain metrics
# restricted to genes that were already identified upstream in the
# pipeline as NS-associated (NS = Nsyn = nonsynonymous-variant enrichment,
# based on the ratio of nonsynonymous to synonymous SNVs/substitutions):
#   - "ns_increased_genes" -> genes with increased nonsynonymous enrichment
#                              (pN/pS or dN/dS) in IBD vs Control
#   - "ns_decreased_genes" -> genes with decreased nonsynonymous enrichment
#                              in IBD vs Control
# For these genes, it computes:
#   1) mean gene count retained per genome before/after the breadth_minCov
#      QC filter (a QC/reporting check, output currently commented out),
#   2) count-weighted mean pN/pS and dN/dS per gene x Phenotype (output
#      currently commented out),
#   3) the correlation between pN/pS and clinical inflammation severity
#      (fecal calprotectin), separately for NS-increased vs NS-decreased
#      genes and by Phenotype -- this is the script's main saved output.
#
# This script is run once per species as a SLURM array job: "file" is the
# path to that species' per-gene inStrain output (.tsv).
#
# EXPECTED OBJECTS IN ENVIRONMENT (set by the calling/wrapper script, not
# defined here):
#   file              - path to the species' gene-level inStrain .tsv
#   genome_presence   - data frame containing qc'ed species-sample pairs
#   gp_path           - path to the genome_presence file
#   samples_phenotype - sample metadata (label -> subject, Phenotype, cohort)
#   ns_increased_genes_list_file - list of genes selected from linear regression model
#   ns_decreased_genes_list_file - list of genes selected from linear regression model
#   output_wd, task   - output directory / sub-folder for this analysis
#
# OUTPUT (written to output_wd/task/):
#   <species>_pnps_corr_severity.csv - Spearman correlation between pN/pS
#       and fecal calprotectin (inflammation severity), per Phenotype and
#       NS gene category (Nsyn-increased / Nsyn-decreased)
#
# (Two additional outputs -- mean gene counts before/after filtering, and
# mean pN/pS & dN/dS per gene/Phenotype -- are computed but their
# write_csv() calls are commented out below, so they are not saved by
# default.)
##############################################################################

# Minimum breadth at min-coverage required for a gene to be considered
# reliably genotyped in a sample.
breadth_minCov_TH <- 0.2

# Extract the species name from the input file path.
sp_name <- str_remove(str_split_fixed(file,"/",n=12)[,12],"\\.tsv")


# Genome presence table for this species: which genome was called present
# in which sample.
genome_presence <- read.delim(file=paste0(gp_path,sp_name,"_genome_sample_presence.csv"), header=T, sep=",") %>% dplyr::select(label,genome)


# list of genes of interest
# Gene sets flagged upstream (from an earlier NS/enrichment model, e.g.
# script 10's output) as having significantly increased or decreased
# nonsynonymous (Nsyn) variant enrichment in IBD vs Control, pooled across
# both IBD subtypes (CD and UC).
ns_increased_genes <- read.delim(ns_increased_genes_list_file, sep="\n", header = F) %>% pull()
print(length(ns_increased_genes))

ns_decreased_genes <- read.delim(ns_decreased_genes_list_file, sep="\n", header = F) %>% pull()
print(length(ns_decreased_genes))

# Column names for inStrain's gene_info.tsv output (see script 2 for
# details on each column).
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
  dplyr::select(-filename_scaffold)

#print(head(sp_gene_info.raw))
#print(nrow(sp_gene_info.raw))

# Build the per-subject, per-gene table *before* applying the
# breadth_minCov QC filter (kept separately so the "before filtering" gene
# counts below can be computed): fill in missing gene lengths, attach
# Phenotype/cohort metadata, keep the most complete instance of each gene
# per subject, and restrict to genome-sample pairs with confirmed genome
# presence.
sp_gene_info.before_filt <- sp_gene_info.raw %>%
  # fill in the missing gene lengths (because later we use them as weights)
  group_by(gene) %>%
  fill(gene_length, .direction = "downup") %>% ungroup() %>%
  inner_join(samples_phenotype, by = "label") %>%
  # select the most complete gene instance per subject # breadth_min_cov #
  group_by(subject,gene) %>%
  arrange(desc(breadth_minCov)) %>%
  slice_head(n=1) %>% ungroup() %>% 
  mutate(Phenotype = factor(Phenotype , levels = c("Control", "CD", "UC")),
  cohort = factor(cohort, levels = c("Pediatric","Adult"))) %>%
  # select samples that have passed the 0.5 breadth threshold
  inner_join(
    dplyr::select(genome_presence, c("label","genome")),
  by = c("label","genome"))

# The QC-filtered table: same as above, plus the breadth_minCov threshold
# applied. This is the table used for all downstream analysis in this
# script.
sp_gene_info <-  sp_gene_info.before_filt %>%
  filter(breadth_minCov >= breadth_minCov_TH )

# mean number of genes per genome we retrieved before and after filtering 
# QC/reporting check: how many genes per genome-sample are retained before
# vs after the breadth_minCov filter, averaged across samples per genome.
# Useful for assessing how much data the QC filter removes for this
# species.
n_genes_before <- sp_gene_info.before_filt %>% 
group_by(genome,label) %>%
summarise(n_genes_before_filt = n()) %>%
group_by(genome) %>%
summarise(mean_n_genes_before_filt = mean(n_genes_before_filt, na.rm =T))

n_genes_after <- sp_gene_info %>% 
group_by(genome,label) %>%
summarise(n_genes_after_filt = n()) %>%
group_by(genome) %>%
summarise(mean_n_genes_after_filt = mean(n_genes_after_filt, na.rm =T))

n_genes_stats <- inner_join(
  n_genes_before,
  n_genes_after,
  by = c("genome"))

# Output currently disabled (commented out) -- kept as an available QC
# export if needed.
#output <- paste0(output_wd, "/", task, "/", sp_name, "_Ngene_mean.csv")
#write_csv(n_genes_stats, file = output , quote = "none", num_threads = 2)


# mean dnds/pnps FOR GENES NS-decreased or NS-increased
# For genes in either the NS-increased or NS-decreased sets, compute the
# count-weighted mean pN/pS (weighted by total SNV nonsyn+syn counts) and
# dN/dS (weighted by total SNS nonsyn+syn counts) per gene x Phenotype.
# Weighting by variant counts gives more influence to genes/samples with
# more variant evidence.
sp_pnps <- sp_gene_info %>%
filter(!is.na(dNdS_substitutions) | !is.na(pNpS_variants)) %>%
filter(gene %in% c(ns_decreased_genes,ns_increased_genes)) %>%
group_by(Phenotype,gene) %>%
 reframe(mean_pNpS = weighted.mean(pNpS_variants, w = SNV_N_count + SNV_S_count, na.rm = T),  
         mean_dNdS = weighted.mean(dNdS_substitutions,  w = SNS_N_count + SNS_S_count, na.rm = T)) %>%
 mutate(genome = sp_name) %>%
 dplyr::select(gene,Phenotype,mean_pNpS,mean_dNdS, genome)

output <- paste0(output_wd, "/", task, "/", sp_name, "_pnps_dnds_mean.csv")
write_csv(sp_pnps, file = output , quote = "none", num_threads = 2)

