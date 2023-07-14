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
plot_dir <- c("~/Documents/whi_sca/rna/plots") # Where to save output plots
id_dir <- ("~/Documents/whi_sca/rna/ids")
meta_dir <- c("~/Documents/whi_sca/rna/meta")
results_dir <- c("~/Documents/whi_sca/rna/results")

data_name <- c("whi_topmed_to6_rnaseq_gene_reads.gct.gz") 

covar <- read.csv(file = paste0(meta_dir, "/sct_all_covars_Jul23.csv"))

# RNA seq data ------------------------------------------------------------------
data <- read_omic(name = data_name, wd = wd)

sub <- as.matrix(data[,colnames(data)%in% covar$rnaseq_ids])
# TODO (turn into an option) # Filter data to just include sct genotyped samples, and remove name/descript cols 
sct_ids <- covar[covar$sct == 1, rnaseq_ids] # 144
sct_covar <- covar[covar$sct == 1, ]

sub <- sub[, colnames(sub) %in% sct_ids] # OR
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
pData <- select(pData, plate, sct) # TODO test with and without sct covar

combat <- sva::ComBat_seq(counts=exprs, batch=pData$plate, full_mod=TRUE, group = pData$sct)
row.names(combat) <- dge$genes$Name
dim(combat)

#### Normalise -----------------------------------------------------------------
combat_norm <- vst(combat)
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

write.csv(pheno_exprs, file = paste0(meta_dir, "/sct_all_covars_hbg_Jul23.csv"), row.names = F)
saveRDS(dge_bc, file = "dge_bc_all.rds")
