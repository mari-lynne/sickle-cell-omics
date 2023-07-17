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
library(readODS)
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
pheno <- filter(pheno, sct == "1" & flagged_seq_qc != "Yes")

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

# Update factors
covar$plate <- as.factor(covar$plate)

# Add whi data
whi <-
  read.csv(file = paste0(meta_dir, "/WHI_HTNKidney_2023-02-27.csv")) %>% 
  rename(subject_id = subjectID)

whi_sct <- filter(whi, subject_id %in% pheno$subject_id)
covar <- left_join(covar, whi_sct, by= "subject_id") %>% select(where(not_all_na))

covar <- rename(covar, eGFR = eGFRCKDEpi)

# Form DGE list ---------------------------------------------------------------

dge <- DGEList(counts=counts, samples=covar, genes=gene_info)

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

# Check that counts are normalised/PCA
lcpm2 <- cpm(dge$counts, log=TRUE)
pdf(paste(plot_dir,"count-filter-sct1.pdf",sep ="/"))
par(mfrow=c(1,2))
plot(density(lcpm[,1:100]), main="Pre-filtering", xlab="log CPM", col="seagreen")
plot(density(lcpm2[,1:100]), main="Post-filtering", xlab="log CPM", col="#ef5a3a")
dev.off()

# filtered count data
lcpm2 <- cpm(dge, log = TRUE)

# Pick plate 
group <- as.factor(dge$samples$plate)
col.group = group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set3")
col.group <- as.character(col.group) # Character vector of colours
# Plot MDS, PCA is faster
plotMDS(lcpm2, labels = group, col = col.group, gene.selection = "common")

# TODO filter two outlier samples
# Although it seems like these are big outliers, the scale is pretty small, include PC's and batch in lm should correct for this


# Normalisation----------------------------------------------------------------

dge <- calcNormFactors(dge, method = "TMM")
# dge  <- estimateDisp(dge, robust=TRUE)

# Design Matrix ---------------------------------------------------------------

# missing smoking recode into never for now

dge$samples$smoking <- ifelse(dge$samples$smoking == "Missing", "Never", dge$samples$smoking)
dge$samples$smoking <- ifelse(dge$samples$smoking == "Never", 0, 1)

design <-
  model.matrix(~eGFR + age + plate + bmi_t0 + smoking +
                 PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10 +PC11 +PC12 +PC13 +PC14,
               data = dge$samples)

des <- model.matrix( ~eGFR, data = dge$samples)

View(design)

# Limma voom -------------------------------------------------------------------

v <- voom(dge, des, plot = TRUE)

vfit <- lmFit(v, design)

# Contrast matrix -------------------------------------------------------------

# Define contrasts to make
contr.mat <- makeContrasts(sct = sct1 - sct0, levels = colnames(design))
vfit <- contrasts.fit(vfit, contrasts=contr.mat)
efit <- eBayes(vfit)
plotSA(efit) # Needs improvement

sct_top <- topTable(efit, coef="eGFR", number=50)

sct_top <- filter(sct_top, P.Value <0.005)

write.csv(sct_top, file = paste0(results_dir,"/sct_linear_genes.csv"))


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

# Haemostasis related genes ----------------------------------------------------

# DGE SCA Steady state vs Healthy controls
sca_genes <- read.csv(file = "sca_dge.csv")
sca_genes$Description <- str_trim(sca_genes$Description)
sig_sca <- filter(toptab, gene_name %in% sca_crisis$Description)

# DGE SCA Crisis vs Healthy controls
sca_crisis <- read.csv(file = "sca_crisis_dge.csv")
# trim white space
sca_crisis$Description <- str_trim(sca_crisis$Description)


a <- plot_vol(toptab, lab = 'adj.sig',
              genes = sca_genes$Description,
              title = "SCT (rs334) vs healthy controls \n\n SCA steady state genes",
              FC = 0.25,
              alpha = 0.6,
              colours = c("orange", "#ef5a3a", "#a50000", "#800000"))


b <- plot_vol(toptab, lab = 'adj.sig',
         genes = sca_crisis$Description,
         title = "SCT (rs334) vs healthy controls \n\n SCA crisis genes",
         FC = 0.25,
         alpha = 0.6,
         colours = c("orange", "#ef5a3a", "#a50000", "#800000"))

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
