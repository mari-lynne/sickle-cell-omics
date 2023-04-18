### Aims:

#1)  Pre-process RNAseq data for DGE analysis of sickle cell trait
#2)  Run differential gene expression analysis at both a global level, and individual gene scale

# Steps before (phenotype summamry statistics), steps after (pathway analysis)

## 1) Pre-processing

### R set up and Functions -----------------------------------------------------

# Bioconductor packages
library(BiocManager)
library(edgeR)
library(limma)

# Data viz
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

# Data cleaning
library(readODS)
library(forcats)
library(stringi)
library(stringr)
library(janitor)
library(data.table)
library(dplyr)
library(tidylog)

source("~/scripts/functions.R")

# Options and Directories ------------------------------------------------------

wd <- c("~/Documents/whi_sca/rna") # Where main data is saved
plot_dir <- c("~/Documents/whi_sca/rna/plots") # Where to save output plots
id_dir <- ("~/Documents/whi_sca/rna/ids")
meta_dir <- c("~/Documents/whi_sca/rna/meta")
results_dir <- c("~/Documents/whi_sca/rna/results")

data_name <- c("whi_topmed_to6_rnaseq_gene_reads.gct.gz") 
pheno_name <- c("sct_rnaseq_pheno.txt")

setwd(wd)

# Data Set up ------------------------------------------------------------------

# RNA seq
data <- read_omic(name = data_name, wd = wd)

# Phenotype and meta data (see sct_id_pheno_match.R for more)
pheno <- fread(paste(meta_dir, pheno_name, sep = "/"))

# Filter data to just include sct genotyped samples 
sub <- data[,colnames(data)%in% pheno$rnaseq_ids]

# Make DGE list object
counts <- as.matrix(sub)

# Sum of mapped reads per sample
lsize <- colSums(counts)
# save in pheno 
pheno$lib.size_og <- lsize

gene_info <- data[, 1:2]

# PCA -------------------------------------------------------------------------

# See pca.R
pc_top <- fread(file = paste0(results_dir, "/pc_top.txt"))

covar <- cbind(pheno, pc_top)
covar <- covar %>% select(-V1)

covar <- covar %>%
  mutate(plate = coalesce(covar$pick_plate1, covar$pick_plate_2_repeat)) %>%
  select(-c("pick_plate1", "pick_plate_2_repeat"))

# Update factors
covar$plate <- as.factor(covar$plate)
covar$sct <- as.factor(covar$sct)

# Form DGE list ---------------------------------------------------------------

dge <- DGEList(counts=counts, samples=covar, genes=gene_info, group = covar$sct)

# QC checks -----------------------
lcpm <- cpm(dge$counts, log=TRUE)

L <- mean(dge$samples$lib.size) * 1e-6
M <- median(dge$samples$lib.size) * 1e-6
c(L, M)

# Filter zero counts across all samples
zero_counts <- rowSums(dge$counts==0)==ncol(dge$counts)
dge <- dge[!zero_counts, ,keep.lib.sizes=FALSE]
dim(dge$counts)

# Filter low count
keep.exprs <- filterByExpr(dge)
dge <- dge[keep.exprs, ,keep.lib.sizes=FALSE]
dim(dge) # 19254   837

lcpm2 <- cpm(dge$counts, log=TRUE)

pdf(paste(plot_dir,"count-filter.pdf",sep ="/"))
par(mfrow=c(1,2))
plot(density(lcpm[,1:100]), main="Pre-filtering", xlab="log CPM", col="seagreen")
plot(density(lcpm2[,1:100]), main="Post-filtering", xlab="log CPM", col="#ef5a3a")

dev.off()

# Check that counts are normalised/PCA

# filtered count data
lcpm2 <- cpm(dge, log = TRUE)
# Set up plate as group/factor to plot
group <- as.factor(dge$samples$plate)
# Make vector of colours for plotting
col.group = group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set3")
col.group <- as.character(col.group) # Character vector of colours
# Plot MDS, PCA is faster
plotMDS(lcpm2, labels = group, col = col.group, gene.selection = "common")
# These could be the flagged QC samples chec
group = as.factor(dge$samples$flagged_seq_qc)
col.group <- group
levels(col.group) <- brewer.pal(nlevels(col.group), "Paired")
col.group <- as.character(col.group)
plotMDS(lcpm2, labels = group, col = col.group, gene.selection = "common")
# Doesn't really look like it

# Pick plate 2 PCA
group <- as.factor(dge$samples$pick_plate_2_repeat)
col.group = group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set3")
col.group <- as.character(col.group) # Character vector of colours
# Plot MDS, PCA is faster
plotMDS(lcpm2, labels = group, col = col.group, gene.selection = "common")
# Not this either

group <- as.factor(dge$samples$sct)
col.group = group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set2")
col.group <- as.character(col.group) # Character vector of colours
# Plot MDS, PCA is faster
plotMDS(lcpm2, labels = group, col = col.group, gene.selection = "common")

# Although it seems like these are big outliers, the scale is pretty small, include PC's and batch in lm should correct for this


# Normalisation----------------------------------------------------------------

dge <- calcNormFactors(dge, method = "TMM")
# dge  <- estimateDisp(dge, robust=TRUE)

lcpm_norm <- cpm(dge, log=TRUE)
boxplot(lcpm_norm, las=2, col=col, main="")
title(main="B. Example: Normalised data", ylab="Log-cpm")

# Design Matrix ---------------------------------------------------------------
# Make var of interest a factor
design <-
  model.matrix(~0 +
                 sct + age + plate + bmi_t0 +
                 PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10 +PC11 +PC12 +PC13 +PC14,
               data = dge$samples)

View(design)

# Limma voom -------------------------------------------------------------------

pdf(paste(plot_dir,"voom2.pdf",sep ="/"))
v <- voom(dge, design, plot = TRUE)
dev.off()

group <- as.factor(dge$samples$plate)
col.group = group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set3")

norm.expr <- v$E
boxplot(norm.expr, 
        col=col.group,
        main="Distribution of normalised counts",
        xlab="",
        ylab="log2 normalised expression",
        las=2,cex.axis=0.8)

vfit <- lmFit(v, design)

# Contrast matrix -------------------------------------------------------------

# Define contrasts to make
contr.mat <- makeContrasts(sct = sct1 - sct0, levels = colnames(design))
vfit <- contrasts.fit(vfit, contrasts=contr.mat)
efit <- eBayes(vfit)
plotSA(efit)

sct_top <- topTable(efit, coef = 1, number = Inf)


# Volcano Plots ---------------------------------------------------------------
# Caused by error in `FUN()`: ! object 'gene_name' not found
# update colname to match function
#TODO improve function to have a gene col vector name variable

toptab <- sct_top
toptab$gene_name <- toptab$Description

# TODO 
# sig function also include FC threshold, check gray shading bit for log FC

pdf(paste(plot_dir,"sct_dge2.pdf",sep ="/"))
plot_vol(toptab, lab = 'top', title = "DEG Sickle Cell Trait (rs334)", FC = 0.25, alpha = 0.6,
          colours = c("orange", "#ef5a3a", "#a50000", "#800000"))
dev.off()

png(paste(plot_dir,"sct_dge2.png",sep ="/"))
plot_vol(toptab, lab = 'top', title = "DEG Sickle Cell Trait (rs334)", FC = 0.25, alpha = 0.6,
         colours = c("orange", "#ef5a3a", "#a50000", "#800000"))
dev.off()

save.image(file = paste0(results_dir, "/sct_dge.RData"))

# sig FC gene list
# Adj P < 0.05 (FDR), and FC >0.02

sig <- filter(sct_top, adj.P.Val <= 0.05)

results_dir <- paste(wd, "results", sep = "/")
write.table(sig, file = paste(results_dir, "top_sig_dge.txt", sep = "/"), quote = F, sep = "\t", row.names = F, col.names = T)
write.table(sig$Description, file = paste(results_dir, "top_genes.txt", sep = "/"), quote = F, sep = "\t", row.names = F, col.names = F)


# Proteomics -------------------------------------------------------------------
# Cross ref proteomics list of genes, replot

proteins <- fread(paste0(wd, "/meta/sct_proteins.txt"))

# Calculate PCA for RNA seq data
# Get BMI covars

png(paste(plot_dir,"sct_dge_proteins2.png",sep ="/"))
plot_vol(toptab, lab = 'sig', genes = proteins$Protein, title = "DEG Sickle Cell Trait (rs334)", FC = 0.25, alpha = 0.6,
         colours = c("orange", "#ef5a3a", "#a50000", "#800000"))
dev.off()

pdf(paste(plot_dir,"sct_dge_proteins2.pdf",sep ="/"))
plot_vol(toptab, lab = 'sig', genes = proteins$Protein, title = "DEG Sickle Cell Trait (rs334)", FC = 0.25, alpha = 0.6,
         colours = c("orange", "#ef5a3a", "#a50000", "#800000"))
dev.off()

sig_list <- filter(toptab, gene_name %in% proteins$Protein)
write.csv(sig_list, file = paste0(results_dir, "/dge_proteiomics.csv"), row.names = F)
