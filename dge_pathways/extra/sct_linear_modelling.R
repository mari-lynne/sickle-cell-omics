# RNAseq for modelling linear traits

# Pipeline based on https://github.com/nkurniansyah/Olivia

# Set up:

source("~/scripts/r-packages.R")
source("~/scripts/functions.R")
source("~/scripts/color-palettes.R")

# Options and Directories ------------------------------------------------------

wd <- c("~/Documents/whi_sca/rna") # Where main data is saved
plot_dir <- c("~/Documents/whi_sca/rna/plots") # Where to save output plots
id_dir <- ("~/Documents/whi_sca/rna/ids")
meta_dir <- c("~/Documents/whi_sca/rna/meta")
results_dir <- c("~/Documents/whi_sca/rna/results")

data_name <- c("whi_topmed_to6_rnaseq_gene_reads.gct.gz") 
pheno_name <- c("sct_rnaseq_pheno.txt")

setwd(wd)

# Load data --------------------------------------------------------------------  

# RNA seq
data <- read_omic(name = data_name, wd = wd)

# Pheno

# Phenotype and meta data (see sct_id_pheno_match.R for more)
pheno <- fread(paste(meta_dir, pheno_name, sep = "/"))
# pheno <- filter(pheno, sct == "1" & flagged_seq_qc != "Yes")
pheno <- filter(pheno, flagged_seq_qc != "Yes")

# Filter data to just include sct genotyped samples 
sub <- data[,colnames(data)%in% pheno$rnaseq_ids]

# Make DGE list object
counts <- as.matrix(sub)

# Genes
gene_info <- data[, 1:2]

# PCA -------------------------------------------------------------------------

# See pca.R
pc_top <- fread(file = paste0(results_dir, "/pc_top.txt"))
pc_top <- rename(pc_top, rnaseq_ids = V1)

covar <- left_join(pheno, pc_top, by = "rnaseq_ids")

covar <- covar %>%
  mutate(plate = coalesce(covar$pick_plate1, covar$pick_plate_2_repeat)) %>%
  select(-c("pick_plate1", "pick_plate_2_repeat"))

# Recode smoking never missing vs current and past
covar <- covar %>%
  mutate(smoking = ifelse
         (smoking == "Missing", "Never", covar$smoking))

# Update factors
covar$plate <- as.factor(covar$plate)

# Add whi data
whi <-
  read.csv(file = paste0(meta_dir, "/WHI_HTNKidney_2023-02-27.csv")) %>% 
  rename(subject_id = subjectID)

whi_sct <- filter(whi, subject_id %in% pheno$subject_id)
covar <- left_join(covar, whi_sct, by= "subject_id") %>% select(where(not_all_na))

covar <- rename(covar, eGFR = eGFRCKDEpi)

# Summarise pheno
summary_phen<- summarize_phenotypes(pheno = covar,
                                    categorical_variables = c("smoking"),
                                    numeric_variables = c("age","eGFR","eGFRMDRD"),
                                    strata = "smoking")
summary_phen

# check ids match
rnaseq_matrix <- counts
rownames(covar) <- covar$rnaseq_ids

IDs_both <- intersect(rownames(covar), colnames(rnaseq_matrix))
rnaseq_matrix <- rnaseq_matrix[, IDs_both]

# rnaseq matrix rownames should be gene_ids
rownames(rnaseq_matrix) <- gene_info$Description

# define trait of interest to study as an exposure associated with genes
trait <- "eGFR"
covariates_string <- "age + as.factor(plate) + as.factor(smoking)+
bmi_t0 +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10 +PC11 +PC12 +PC13 +PC14"
  
resid_plot<- residual_plot(pheno = covar, 
                           traits = trait,
                           covariates_string = covariates_string)


# QC and plots --------------------------------------------------------------
# clean counts 
median_norm <- median_normalization(rnaseq_matrix)

xrange <- range(colSums(rnaseq_matrix))
par(mfrow = c(1,2))
hist(colSums(rnaseq_matrix), 
     xlim = xrange, 
     main = "",
     xlab = "Before median normalization",
     ylab = "Library size distribution")
hist(colSums(median_norm), 
     xlim = xrange, 
     main = "",
     xlab = "After median normalization",
     ylab = "Library size distribution")

clean_count_matrix <- apply_filters(count_matrix = median_norm, 
                                    median_min = 1, 
                                    expression_sum_min = 10, 
                                    max_min = 10, 
                                    range_min = 5, 
                                    prop_zero_max = 0.5)

par(mfrow=c(1,2))
plot(density(log_replace_half_min(median_norm)[,3]),
     main="Before filtering", xlab="log trascript")
plot(density(log_replace_half_min(clean_count_matrix)[,3]), 
     main="After filtering", xlab="log trascript")

log_counts<- log_replace_half_min(clean_count_matrix)
log_counts<- melt(log_counts)

box_plot<- ggplot(log_counts, aes(x = Var2, y = value)) + 
  stat_boxplot(aes(Var2, value), 
               geom='errorbar', linetype=1, width=0.5)+ 
  xlab("Sample")+ 
  ylab("log(Transcripts)")+ 
  geom_boxplot( aes(Var2, value),outlier.shape=1)+
  stat_summary(fun = mean, geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dotted") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

box_plot


# Log TPM ----------------------------------------------------------------------





# TWAS -------------------------------------------------------------------------

set.seed(12)

# rownames of pheno/covar must match colnames counts
# geneIDs must be colnames of count_matrix or specified

quantile_emp_trascript<-lm_count_mat_emp_pval(clean_count_matrix, 
                                              pheno=covar, 
                                              trait=trait,
                                              covariates_string=covariates_string,
                                              n_permute=80,
                                              log_transform = "log_replace_half_min",
                                              outcome_type ="continuous",
                                              gene_IDs=NULL)

tophits <- quantile_emp_trascript[which(quantile_emp_trascript$bh_emp_pvals< 0.05),]
head(tophits)

# no better including nSCT

## Specific genes --------------------------------------------------------------

dge <- read.csv(file = paste0(results_dir, "/top_dge_genes.csv"))
gene_names <- dge$Description

set <- intestect(gen_names) # dont use gene names not valid i think stick with esn to avoid dup issue when merging

set.seed(12)
perm_res<- lm_count_mat_perm_pval(count_matrix=clean_count_matrix,
                                  pheno=covar,
                                  trait=trait, 
                                  covariates_string=covariates_string,
                                  n_permute=100000,
                                  gene_IDs=gene_names,
                                  seed = NULL,
                                  log_transform = "log_replace_half_min",
                                  outcome_type ="continuous")

# however in general these methods are just measuring individual gene significance
# with effect (beta coefficient) of our linear variable

# we want to see if a set of genes together could be used to predict outcome y
# preprocess data as is input into these linear models 
# except we use them as variables in feature selection/modelling
