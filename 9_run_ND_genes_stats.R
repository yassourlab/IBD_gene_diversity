##############################################################################
# 9_run_ND_genes_stats.R
#
# PURPOSE
# -------
# For a single species, pull out the per-sample, per-gene inStrain metrics
# (nucleotide diversity, dN/dS, pN/pS, SNV/SNS counts) restricted to genes
# that were already identified upstream in the pipeline as ND-associated:
#   - "purged_genes"      -> genes with decreased nucleotide diversity (ND)
#                             in IBD relative to Control ("purged" of
#                             diversity, e.g. under purifying/directional
#                             selection)
#   - "diversified_genes" -> genes with increased nucleotide diversity (ND)
#                             in IBD relative to Control
# This produces a clean, gene-filtered table (one row per gene x sample)
# that downstream scripts/plots use to compare these two ND-defined gene
# categories.
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
#   output_wd, task   - output directory / sub-folder for this analysis
#
# OUTPUT (written to output_wd/task/):
#   <species>_dnds.csv - per-sample, per-gene table (ND, dN/dS, pN/pS, and
#                        SNV/SNS synonymous/nonsynonymous counts) for genes
#                        in the purged or diversified ND gene sets
##############################################################################

# Minimum breadth at min-coverage required for a gene to be considered
# reliably genotyped in a sample.
breadth_minCov_TH <- 0.2

# Extract the species name from the input file path.
sp_name <- str_remove(str_split_fixed(file,"/",n=12)[,12],"\\.tsv")


# Genome presence table for this species (as in earlier scripts): which
# genome was called present in which sample.
genome_presence <- read.delim(file=paste0(gp_path,sp_name,"_genome_sample_presence.csv"), header=T, sep=",") %>% dplyr::select(label,genome)

# list of genes of interest
# Gene sets flagged upstream (from an earlier ND model, e.g. script 2's
# output) as having significantly decreased ("purged") or increased
# ("diversified") nucleotide diversity in IBD vs Control, pooled across
# both IBD subtypes (CD and UC).
purged_genes <- read.delim("/sci/backup/morani/lab/Projects/IBD_diversity/purged_genes_both_IBD.txt", sep="\n", header = F) %>% pull()
print(length(purged_genes))
diversified_genes <- read.delim("/sci/backup/morani/lab/Projects/IBD_diversity/diversified_genes_all.txt", sep="\n", header = F) %>% pull()
print(length(diversified_genes))


# Column names for inStrain's gene_info.tsv output (see script 2 for
# details on each column).
column_names <- c("filename_scaffold", "gene","gene_length",
"coverage","breadth","breadth_minCov","nucl_diversity","start",
"end","direction","partial","dNdS_substitutions","pNpS_variants",
"SNV_count","SNV_S_count","SNV_N_count","SNS_count","SNS_S_count","SNS_N_count",
"divergent_site_count")


# Read and parse the raw per-gene inStrain output for this species (same
# parsing/genome-ID cleanup logic as in scripts 2 and 4).
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
  # select the most complete gene instance per subject # breadth_min_cov #
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

# FOR GENES ND-decreased or ND-increased
# Restrict to genes that have at least one of dN/dS or pN/pS computed, and
# that belong to either the purged (ND-decreased) or diversified
# (ND-increased) gene sets. Keep only the columns needed for downstream
# comparison of these two gene categories.
sp_dnds <- sp_gene_info %>%
filter(!is.na(dNdS_substitutions) | !is.na(pNpS_variants)) %>%
#filter(dNdS_substitutions > 0 | pNpS_variants > 0) %>%
filter(gene %in% c(purged_genes,diversified_genes)) %>%
dplyr::select(label,gene,genome,gene_length,breadth_minCov,nucl_diversity,dNdS_substitutions,pNpS_variants,SNV_S_count,SNV_N_count,SNS_S_count,SNS_N_count)
#group_by(label) %>%


print(nrow(sp_dnds))

output <- paste0(output_wd, "/", task, "/", sp_name, "_dnds.csv")
write_csv(sp_dnds, file = output , quote = "none", num_threads = 2)

