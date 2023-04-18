### R set up and Functions -----------------------------------------------------

# Bioconductor packages
library(BiocManager)
library(edgeR)
library(limma)

library(biomaRt)
library(Organism.dplyr)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(DOSE)
library(org.Hs.eg.db)

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

# Pathway analysis
library(clusterProfiler)

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

load(file = paste0(results_dir, "/sct_dge.RData"))

table(pheno$age)

# 1) Get list of sig genes -----------------------------------------------------

top_dge <- toptab %>% filter(adj.P.Val <= 0.05 & (logFC < -0.2 | logFC > 0.2))
dge_genes <- unique(top_dge$gene_name)
up_dge <- toptab %>% filter(adj.P.Val <= 0.05 & logFC > 0.2) # 29
down_dge <- toptab %>% filter(adj.P.Val <= 0.05 & logFC < -0.2) # 29

write.csv(top_dge, file = paste0(results_dir, "/top_dge_genes.csv"))
write.csv(toptab, file = paste0(results_dir, "/toptab.csv"))

# Format for clusterprofiler R

# Run Pathway Enrichment analysis ----------------------------------------------

# 2) Set up gene universe ------------------------------------------------------

src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")

# Search keys
chrom <- c(as.character(1:22))
chrom <- str_c(rep("chr", length(chrom)), chrom)
genes <- dge_genes

universe <- Organism.dplyr::select(src, 
                   keys = chrom,
                   columns = c("symbol","entrez"),
                   keytype = "cds_chrom")

universe <- distinct(universe)

genes_ez <- Organism.dplyr::select(src, 
                   keys = genes,
                   columns = c("symbol","entrez"),
                   keytype = "symbol")

genes_ez <- distinct(genes_ez)

# 3) run GO code ---------------------------------------------------------
ego <- enrichGO(gene = genes_ez$entrez,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", readable = TRUE)
View(ego@result)


edo <- enrichDO(gene = genes_ez$entrez,
                ont = "DO",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", readable = TRUE)

View(edo@result)

go_dge <-
  groupGO(
    gene = dge_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    level = 3,
    readable = TRUE
  )


# MsigDB

c7 <- read.gmt("~/Downloads/public_data/c7.all.v2023.1.Hs.entrez.gmt")
c5 <- read.gmt("~/Downloads/public_data/c5.all.v2023.1.Hs.entrez.gmt")

c7_ego <- enricher(genes_ez$entrez, TERM2GENE=c7) #c5
c7_sig <- subset(c7_ego@result, p.adjust <0.05)

View(c7_ego@result)

dotplot(c7_ego, showCategory=15)



# Data viz ---------------------------------------------------------------------

