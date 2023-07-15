# 1_SCT RNA seq pre-process:

# Aims:
# Gather rnaseq counts, and covariate data -> dge list
# Check for sample outliers
# Combat adjustment for plate #TODO visulasation of combat
# Normalise
# Calculate HBB ratios
# Calculate transcripts per million (TPM)
# 
# Split data for both just SCT and all samples

# Set up -----------------------------------------------------------------------

source("~/scripts/r-packages.R")
source("~/scripts/functions.R")

wd <- c("~/Documents/whi_sca/rna") # Where main data is saved
plot_dir <- c("~/Documents/whi_sca/rna/plots")
meta_dir <- c("~/Documents/whi_sca/rna/meta")
results_dir <- c("~/Documents/whi_sca/rna/results")

data_name <- c("whi_topmed_to6_rnaseq_gene_reads.gct.gz") 
covar <- read.csv(file = paste0(meta_dir, "/sct_all_covars_Jul23.csv"))

# RNA seq data ------------------------------------------------------------------
data <- read_omic(name = data_name, wd = wd)

sub <- as.matrix(data[,colnames(data)%in% covar$rnaseq_ids])
dim(sub)
# OPTION:
# Filter out hispanics
# Filter data to just include sct genotyped samples, and remove name/descript cols
sct_ids <- covar[covar$ethnic==3, "rnaseq_ids"]
sct_ids <- covar[covar$sct == 1, "rnaseq_ids"] # 144

sub <- sub[, colnames(sub) %in% sct_ids] # OR
dim(sub)
sub <- sub # ALL SAMPLES

## Genes and counts ------------------------------------------------------------
# Genes
gene_info <- data[, 1:2]

# More gene info
gene_info <- rename(gene_info, hgnc_symbol = Description)

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

genes <- getBM(
  attributes=c('ensembl_gene_id', 'hgnc_symbol','start_position','end_position'),
  mart = ensembl)

genes$size = genes$end_position - genes$start_position

# Merge ensembl info and current gene info
# TODO update to filter duplicates by lowest expressing transcript
filter_genes <- 
  left_join(gene_info %>% group_by(hgnc_symbol) %>% mutate(id = row_number()),
            genes %>% group_by(hgnc_symbol) %>% mutate(id = row_number()), 
            by = c("hgnc_symbol", "id"))

# Make DGE list
row.names(sub) <- filter_genes$Name
covar <- covar[covar$rnaseq_ids %in% sct_ids,]
dge <- DGEList(counts=sub, samples=covar, genes=filter_genes) # covar for all samps

# Filter genes missing ensembl and transcript length info 
dge <- dge[!is.na(dge$genes$ensembl_gene_id), ]
dim(dge)

# QC --------------------------------------------------------------------------

## Filter low counts

# Filter zero counts across all samples
zero_counts <- rowSums(dge$counts==0)==ncol(dge$counts)
dge <- dge[!zero_counts, ,keep.lib.sizes=FALSE]
dim(dge$counts)

table(dge$samples$sct)
# Filter low count
keep.exprs <- filterByExpr(dge, group = dge$samples$sct)
dge <- dge[keep.exprs, ,keep.lib.sizes=FALSE]
dim(dge) # 13938  144

## Batch correction ------------------------------------------------------------

# TODO add PCA plots
exprs <- as.matrix(dge$counts)
pData <- dge$samples
pData$sct <- as.factor(pData$sct)
pData <- select(pData, plate, sct) 
combat <- sva::ComBat_seq(counts=exprs, batch=pData$plate, full_mod=TRUE, group = pData$sct)

row.names(combat) <- dge$genes$Name
dim(combat)

## Combat PCA plots -----------------------------------------------

# Check technical and biological factors

# Calculate PCA before batch correction
pca_before <- prcomp(t(exprs))
pca_df_before <- as.data.frame(pca_before$x[, 1:2]) %>%
  mutate(Batch = as.factor(pData$plate), sct = as.factor(pData$sct))

# Calculate PCA after batch correction
pca_after <- prcomp(t(combat))
pca_df_after <- as.data.frame(pca_after$x[, 1:2]) %>%
  mutate(Batch = as.factor(pData$plate), sct = as.factor(pData$sct))

# Plot PCA before batch correction
a <- ggplot(pca_df_before, aes(x = PC1, y = PC2, color = sct)) +
  geom_point() +
  stat_ellipse(aes(fill = sct), alpha = 0.6) +
  labs(title = "PCA Before Batch Correction")
a

# Plot PCA after batch correction
b <- ggplot(pca_df_after, aes(x = PC1, y = PC2, color = sct)) +
  geom_point() + 
  stat_ellipse(aes(fill = sct), alpha = 0.6) +
  labs(title = "After Batch Correction")


library(patchwork)
(a+b) + plot_layout(guides = "collect")
# I think plate might be more useful

ggsave("plate_bc_pca.png")

#### Normalise -----------------------------------------------------------------

# OPTION: # Normalise og data (save as dge_bc this round to reuse code)
combat_norm <- vst(combat)
combat_norm <- vst(dge$counts)
# Update DGE object
dge_bc <- DGEList(counts=combat_norm, samples=dge$samples, genes=dge$genes)


## HBF Phenotype ---------------------------------------------------------------

# Subset data for just Hemoglobin expression data
idx <- (dge_bc$genes$hgnc_symbol %in% c("HBB", "HBG1", "HBG2")) 
sub <- dge_bc[idx,]
dim(sub)

# Make table with HBB and HBG counts, ratios and SCT status for plotting
pheno <- sub$samples
genes <- t(sub$genes)
counts2 <- as.data.frame(t(sub$counts))
colnames(counts2) <- genes[2,]

# Calculate ratios
# For now just aggregate HBG1/2
counts2$HBG <- counts2$HBG1 + counts2$HBG2
counts2$Ratio_HBG <- counts2$HBB / counts2$HBG
pheno_exprs <- cbind(dge_bc$samples, counts2)
dim(pheno_exprs)

pheno_exprs <- clean_cols(pheno_exprs)
pheno_exprs <- select(pheno_exprs, -c("baa23", "as311", "as315", "baa23.1", "as311.1", "as315.1", "whills"))

write.csv(pheno_exprs, file = paste0(meta_dir, "/sct_black_covars_hbg_Jul23.csv"), row.names = F)

# OPTION
# update DGE
dge_bc$samples <- pheno_exprs
saveRDS(dge_bc, file = "dge_bc_black.rds")

# No batch effect DGE
saveRDS(dge_bc, file = "dge_black.rds")
