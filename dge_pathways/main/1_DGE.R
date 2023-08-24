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
library(patchwork)
library(ggpubr)

# Data cleaning
library(forcats)
library(stringi)
library(stringr)
library(janitor)
library(data.table)
library(dplyr)
library(tidylog)

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

# Data Set up ------------------------------------------------------------------
# load(file = paste0(results_dir, "/sct_dge.RData"))

# RNA seq
data <- read_omic(name = data_name, wd = wd)

# Phenotype and meta data (see sct_id_pheno_match.R for more)
pheno <- fread(paste(meta_dir, pheno_name, sep = "/"))

# Filter data to just include sct genotyped samples
pheno <- pheno[flagged_seq_qc != "Yes", ]
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
pc_top <- rename(pc_top, rnaseq_ids = V1)
covar <- left_join(pheno, pc_top, by = "rnaseq_ids")

# Update factors
covar <- covar %>%
  mutate(plate = coalesce(covar$pick_plate1, covar$pick_plate_2_repeat)) %>%
  select(-c("pick_plate1", "pick_plate_2_repeat"))
covar$plate <- as.factor(covar$plate)
covar$sct <- as.factor(covar$sct)

# Updated covar file ----------------------------------------------------------
meta_dir <- c("~/Documents/whi_sca/rna/meta")
samp_covar <- read.csv(file.path(meta_dir, "sct_all_covars_Jul23.csv"))
ciber <- read.csv(file = "~/Documents/whi_sca/rna/results/ciber/R/pp_hat_ciber_merged.csv")
ciber <- rename(ciber, rnaseq_ids = X)

# Include smokestaus and cell fractions
samp_covar <- select(samp_covar, rnaseq_ids, smokestatus)

# Join covar tables :)
samp_covar <- left_join(samp_covar, ciber, by = "rnaseq_ids")
covar <- left_join(covar, samp_covar, by = "rnaseq_ids")

# Form DGE list ---------------------------------------------------------------

dge <- DGEList(counts=counts, samples=covar, genes=gene_info, group = covar$sct)

# QC checks -----------------------
lcpm <- cpm(dge$counts, log=TRUE)
L <- mean(dge$samples$lib.size) * 1e-6
M <- median(dge$samples$lib.size) * 1e-6
c(L, M)

# Filter zero counts across all samples
dim(dge$counts)
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

group <- as.factor(dge$samples$sct)
col.group = group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set2")
col.group <- as.character(col.group) # Character vector of colours
# Plot MDS, PCA is faster
plotMDS(lcpm2, labels = group, col = col.group, gene.selection = "common")

# Although it seems like these are big outliers, the scale is pretty small, include PC's and batch in lm should correct for this

# Outlier removal (see PCA.R) ------------------------------------------------


# Normalisation----------------------------------------------------------------

dge <- calcNormFactors(dge, method = "TMM")
# dge  <- estimateDisp(dge, robust=TRUE)
lcpm_norm <- cpm(dge, log=TRUE)
boxplot(lcpm_norm, las=2, col=col, main="")
title(main="B. Example: Normalised data", ylab="Log-cpm")

# Design Matrix ---------------------------------------------------------------
# Make var of interest a factor
 # just have current versus past/na
# Assuming 'dge' is your data frame
dge$samples$smokestatus <- dge$samples$smokestatus.y
dge$samples$smokestatus <- ifelse(is.na(dge$samples$smokestatus), 0, dge$samples$smokestatus)
dge$samples <- dplyr::mutate(dge$samples, smokestatus = ifelse(smokestatus != 2, 0, smokestatus))
dge$samples$smokestatus <- as.factor(dge$samples$smokestatus)

design <-
  model.matrix(~0 +
                 sct + plate + smokestatus + bmi_t0 + age + 
                 PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10 +PC11 +PC12 +PC13 +PC14 + neut + B_cells + mono,
               data = dge$samples)

# Limma voom -------------------------------------------------------------------

pdf(paste(plot_dir,"voom.pdf",sep ="/"))
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
toptab <- sct_top
toptab$gene_name <- toptab$Description

genes <- c("HEMGN", "RBM38", "SCOC", "ABCB6", "FHDC1", "SPTA1", "TAL1", "RANBP10", "BLVRB", "SOX6", "TRIM10", "YPEL4", "IGF2BP2", "RUNDC3A", "PRDX6", "TRAT1", "SCOC", "BLVRB", "VWCE", "RPL30", "RPS15A")

# RPL30 RPS15A dimond blackfan anemia

pdf(paste(plot_dir,"sct_dge.pdf",sep ="/"))

plot_vol(toptab,
         lab = 'adj.sig',
         title = "DEG Sickle Cell Trait (rs334)",
         FC = 0.1,
         alpha = 0.6,
         gene_lab = "custom",
         genes_interest = genes,
          colours = c("orange", "#ef5a3a", "#a50000", "#800000"))

dev.off()

png(paste(plot_dir,"sct_dge2.png",sep ="/"))
plot_vol(toptab, lab = 'top', title = "DEG Sickle Cell Trait (rs334)", FC = 0.25, alpha = 0.6,
         colours = c("orange", "#ef5a3a", "#a50000", "#800000"))
dev.off()

# save.image(file = paste0(results_dir, "/sct_dge.RData"))

# sig FC gene list
# Adj P < 0.05 (FDR), and FC >0.02

sig <- filter(sct_top, adj.P.Val <= 0.05 & abs(logFC) > 0.15)
sig <- filter(sct_top, adj.P.Val <= 0.05)

write.csv(sig,
          file = paste(results_dir, "top_sig_dge_fc0.csv",
                  sep = "/"),
          row.names = F)

# save.image(file = "sct_july_dge.RData")

write.table(sig$Description, file = paste(results_dir, "top_genes.txt", sep = "/"), quote = F, sep = "\t", row.names = F, col.names = F)


# TOPMED genes of interest -----------------------------------------------------

topmed <- dge[dge$genes$Description %in% c("RUNDC3A", "RUNDC3A-AS1", "HBA1", "HBA2", "HBB", "HBG1", "HBG2", "SLC4A1", "SLC25A39"),]

topmed_samp <- as.data.frame(topmed$samples)
topmed_counts <- as.data.frame(t(topmed$counts))
topmed <- cbind(topmed_samp, topmed_counts)
#topmed <- dplyr::rename(topmed, "TAL1_counts" = "1502")

colnames(topmed)[31:39] <- c("RUNDC3A", "RUNDC3A-AS1", "HBA1", "HBA2", "HBB", "HBG1", "HBG2", "SLC4A1", "SLC25A39")
topmed <- select(topmed, -commonid)

write.csv(topmed, file = paste0(results_dir, "genes-interest.csv"), row.names = F)


# Haemostasis related genes ----------------------------------------------------

# DGE SCA Steady state vs Healthy controls
sig$gene_name <- sig$Description

sca_genes <- read.csv(file = "sca_dge.csv")
sca_genes$Description <- str_trim(sca_genes$Description)
sig_sca <- filter(sig, gene_name %in% sca_genes$Description)

# DGE SCA Crisis vs Healthy controls
sca_crisis <- read.csv(file = "sca_crisis_dge.csv")
sca_crisis$Description <- str_trim(sca_crisis$Description)
sig_crisis <- filter(sig, gene_name %in% sca_crisis$Description)


a <- plot_vol(toptab, lab = 'adj.sig',
              genes = sca_genes$Description,
              gene_lab = "custom",
              title = "SCT (rs334) vs healthy controls \n\n SCA steady state genes",
              FC = 0.25,
              alpha = 0.6,
              colours = c("orange", "#ef5a3a", "#a50000", "#800000"))

a

b <- plot_vol(toptab,
         genes = sca_crisis$Description,
         lab_type = "custom",
         title = "SCT (rs334) vs healthy controls \n\n SCA crisis genes",
         FC = 0.1,
         alpha = 0.6,
         colours = c("orange", "#ef5a3a", "#a50000", "#800000"))

b

(a|b) + plot_annotation(tag_levels = "a")

?patchwork

sig_sca <- filter(toptab, gene_name %in% sca_crisis$Description)


png(paste(plot_dir,"sct_sca_dge_overlap.png", sep ="/"), width = 680, height = 680)
(a|b) 
dev.off()

pdf(paste(plot_dir,"sct_sca_dge_overlap.pdf", sep ="/"))
(a|b)
dev.off()

?png

# TAL1 -------------------------------------------------------------------------

tal1 <- dge[dge$genes$Description == "TAL1",]


tal1_samp <- as.data.frame(tal1$samples)
tal1_counts <- as.data.frame(t(tal1$counts))
tal1 <- cbind(tal1_samp, tal1_counts)
tal1 <- dplyr::rename(tal1, "TAL1_counts" = "1502")

autumn2 <- c("#A92D2F","#CE7726")
comp <- list(c("0", "1"))

tal1$cpm <- tal1$TAL1_counts/1e6

c <-
  tal1 %>% ggplot(aes(x = sct, y = TAL1_counts)) +
  geom_boxplot(fill = autumn2)+
  theme_light() +
  scale_y_continuous(trans='log10') +
  labs(y = "TAL1 expression (counts)\n", x = "SCT status") +
  stat_compare_means(comparisons = comp)

png(paste(plot_dir,"TAL1_exprs.png", sep ="/"))
c
dev.off()

# Methlyation ------------------------------------------------------------------

sub <- dge[dge$genes$Description == "HBB",]

sub_samp <- as.data.frame(sub$samples)
sub_counts <- as.data.frame(t(sub$counts))
sub <- cbind(sub_samp, sub_counts)
sub$sub_counts <- sub[,31]
sub$cpm <- sub$sub_counts/1e6

c <-
  sub %>% ggplot(aes(x = sct, y = sub_counts)) +
  geom_boxplot(fill = autumn2)+
  theme_light() +
  scale_y_continuous(trans='log10') +
  labs(y = "HBB expression (counts)\n", x = "SCT status") +
  stat_compare_means(comparisons = comp)
c

png(paste(plot_dir,"HBB_exprs.png", sep ="/"))
c
dev.off()

# All sig meth genes overlap vol

meth <- read.csv(file = "sct_meth_genes.csv")
meth$Gene <- str_trim(meth$Gene)

m <- plot_vol(toptab, lab = 'sig',
              genes = meth$Gene,
              title = "SCT (rs334) vs healthy controls \n\n Methylation sites",
              FC = 0.25,
              alpha = 0.6,
              colours = c("orange", "#ef5a3a", "#a50000", "#800000"))
m

# Proteomics -------------------------------------------------------------------
# Cross ref proteomics list of genes, replot

proteins <- fread(paste0(wd, "/meta/sct_proteins.txt"))
proteins$Protein <- str_trim(proteins$Protein)

p <- plot_vol(toptab, lab = 'sig',
              genes = proteins$Protein,
              title = "SCT (rs334) vs healthy controls \n\n Proteomics",
              FC = 0.25,
              alpha = 0.6,
              colours = c("orange", "#ef5a3a", "#a50000", "#800000"))
p

pdf(paste(plot_dir,"sct_dge_proteins2.pdf",sep ="/"))
p
dev.off()

png(paste(plot_dir,"sct_dge_proteins.png",sep ="/"))
p
dev.off()

sig_list <- filter(toptab, gene_name %in% proteins$Protein)
write.csv(sig_list, file = paste0(results_dir, "/dge_proteiomics.csv"), row.names = F)

# Protein/meth combo plot

png(paste(plot_dir,"sct_protein_meth.png",sep ="/"))
(p|m)
dev.off()

pdf(paste(plot_dir,"sct_protein_meth.pdf",sep ="/"), paper = "a4r")
(p|m)
dev.off()

## HB- beta to gamma ratios ----------------------------------------------------

# Subset data for just Hemoglobin expression data
HB <- filter(toptab, str_starts(toptab$Description, "HB"))
idx <- (dge$genes$Description %in% c("HBB", "HBG1", "HBG2")) 
sub <- dge[idx,]

# Make table with HBB and HBG counts, ratios and SCT status for plotting

pheno <- sub$samples
genes <- t(sub$genes)
counts2 <- as.data.frame(t(sub$counts))
colnames(counts2) <- genes[2,]

# Calculate ratios
counts2$Ratio_HB_G1 <- counts2$HBB / counts2$HBG1
counts2$Ratio_HB_G2 <- counts2$HBB / counts2$HBG2

# Aggregate HBG (only a base diff so according to Vijay, combine)
# Should I normalize and then add? check scales

a <- pheno_exprs %>%
  ggplot(aes(y= HBG1)) + geom_histogram(bins = 50) +
  ylim(0, 80000) + xlim(0,125) +
  ylab("HBG1 (counts)") + xlab("count (n)") +
  theme_light()

b <- pheno_exprs %>%
  ggplot(aes(y= HBG2)) + geom_histogram(bins = 50) +
  ylim(0, 80000) + xlim(0,125) +
  ylab("HBG2 (counts)") + xlab("count (n)") +
  theme_light()

ggsave(paste(plot_dir,"HBG_hist.png", sep ="/"), width = 3.25, height = 3.25, dpi = 300)
# TODO reinstall ragg - to fix: Error in f(...) : Graphics API version mismatch

# save plot
png(paste(plot_dir,"HBG_hist.png", sep ="/"), width = 1680, height = 1680, res = 300)
(a|b) + plot_annotation(tag_levels = 'A')
dev.off()

# For now just add
counts2$HBG <- counts2$HBG1 + counts2$HBG2
counts2$Ratio_HBG <- counts2$HBB / counts2$HBG
pheno_exprs <- cbind(pheno, counts2)

# Plot ratios :)

comp <- list(c("0", "1"))

options(scipen=10000)
p <- 
  pheno_exprs %>% ggplot(aes(x = group, y = Ratio_HBG, fill = group)) + geom_boxplot() +
  stat_compare_means(comparisons = comp, paired = FALSE, method = "wilcox.test") +
  ylab("HBB/HBG \n") + scale_y_log10() + xlab("SCT") +
  scale_fill_manual(values = c("#A72020", "#FEA82F")) +
  theme_light() + theme(legend.position = "none")


png(paste(plot_dir,"HBG_ratio.png", sep ="/"), width = 1680, height = 1680, res = 300)
p
dev.off()

# HBB ratios voom ----------------------------
# Try using E values from voom (log2cpm)CPM based on normalised library size
idx <- (v$genes$Description %in% c("HBB", "HBG1", "HBG2"))
sub <- v[idx,]

pheno <- sub$targets
genes <- t(sub$genes)
counts2 <- as.data.frame(t(sub$E))
colnames(counts2) <- genes[2,]

# Calc ratios
counts2$HBG <- counts2$HBG1 + counts2$HBG2
counts2$Ratio_HBG <- counts2$HBB / counts2$HBG
pheno_exprs <- cbind(pheno, counts2)

# Plot
comp <- list(c("0", "1"))

options(scipen=10000)
p <-
  pheno_exprs %>% ggplot(aes(x = group, y = Ratio_HBG, fill = group)) + geom_boxplot() +
  stat_compare_means(comparisons = comp, paired = FALSE, method = "wilcox.test") +
  ylab("HBB/HBG \n") + xlab("SCT") +
  scale_y_continuous(trans='log2') +
  scale_fill_manual(values = c("#A72020", "#FEA82F")) +
  theme_light() + theme(legend.position = "none")

png(paste(plot_dir,"HBG_ratio_log2cpm.png", sep ="/"), width = 1680, height = 1680, res = 300)
p
dev.off()


### TPM ---------------------------------------------------------------

# 2) Calculate TPM matrix ---------------------------------------------------------
# function
tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

HB <- c("HBA1", "HBA2", "HBB", "HBG1", "HBG2")
# sub_list = c("TAL1","RUNDC3A", "RUNDC3A-AS1", "HBA1", "HBA2", "HBB", "HBG1", "HBG2", "SLC4A1", "SLC25A39")
sub_list <- filter(toptab, (adj.P.Val < 0.05 & (logFC < -0.1 | logFC > 0.1 ))| Description %in% HB) %>%
  select(Description) # 299 - 162

sub_list <- select(toptab, Description)

sub_list$hgnc_symbol <- sub_list$Description
# Get gene/transcript length from biomart
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(
  attributes=c('ensembl_gene_id', 'hgnc_symbol','start_position','end_position'),
  mart = ensembl)
genes$size = genes$end_position - genes$start_position

sub_list <- inner_join(sub_list, genes, by = "hgnc_symbol")
# filter dups
sub_list <- sub_list[(stri_duplicated(sub_list$hgnc_symbol)==FALSE), ]
dim(sub_list)
topmed <- dge[dge$genes$Description %in% sub_list$Description, ]
dim(topmed)
topmed$genes <- left_join(sub_list, topmed$genes, by = "Description")

  
options(scipen=999)
top_tpm <- as.data.frame(t(tpm(topmed$counts, topmed$genes$size)))
colnames(top_tpm) <- topmed$genes$hgnc_symbol
top_tpm <- cbind(topmed$samples, top_tpm)
top_tpm <- select(top_tpm, -commonid)

# load(file = paste0(results_dir, "/sct_dge.RData"))
# write.csv(toptab, file = paste0(results_dir, "/top_tab.csv"), row.names = F)
# write.csv(top_toptab, file = paste0(results_dir, "/top_tab_sig.csv"), row.names = F)
# write.csv(top_tpm, file = paste0(results_dir, "/tpm_sig.csv"), row.names = F)


# TPM correlations ------------------------------------------------------------
options(scipen=999)
top_tpm <- as.data.frame(t(tpm(topmed$counts, topmed$genes$size)))
colnames(top_tpm) <- topmed$genes$hgnc_symbol #sub_list$
control_vec <- top_tpm$HBB
top_tpm <- as.matrix(top_tpm)
# control_vec <- counts2$Ratio_HBG # This is a proxy for feotal haemoglobin

# Calculate correlation
correlations <- apply(top_tpm, 2, function(x) cor(x, control_vec))
# Output correlations
print(correlations)

# Calculate sig correlation
correlations <- apply(top_tpm, 2, function(x) {
  cor_result <- cor(x, control_vec, method = "spearman")
  
  # Filter correlations based on threshold and statistical significance
  if (abs(cor_result) > 0.5 && cor_result >= 0 && cor_result < 1) {
    p_value <- cor.test(x, control_vec, method = "spearman")$p.value
    if (p_value < 0.05) {
      return(cor_result)
    }
  }
  
  return(NA)  # If correlation doesn't meet the criteria, return NA
})

# Remove NA values
cor_result <- correlations[!is.na(correlations)]

# Output correlations
print(correlations)
