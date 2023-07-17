# ClusterProfiler

# Aims: 
# Set up data for gene enrichment analysis

# Set up -----------------------------------------------------------------------

# Bioconductor packages
library(BiocManager)
library(biomaRt)
library(Organism.dplyr)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(DOSE)

# Data viz
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(patchwork)

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
library(clusterProfiler) # https://yulab-smu.top/biomedical-knowledge-mining-book/index.html
library(org.Hs.eg.db)
library(enrichplot)

# Set up:
source("~/scripts/r-packages.R")
source("~/scripts/functions.R")

## Directories ------------------------------------------------------

wd <- "~/Documents/whi_sca/rna" # Where main data is saved
results_dir <- "~/Documents/whi_sca/rna/results/wgcna"
plot_dir <- "~/Documents/whi_sca/rna/plots/wgcna/updates/just_sct/no_bc" # Where to save output plots
path_dir <- "~/Downloads/public_data/c5" # Download from msigDB
setwd(wd)

## Input ------------------------------------------------------
# Input data should be a list of gene symbols or ensembl/enterez ids to match against db
name <- "/gene_list_violet.csv"
data <- read.csv(file = paste0(results_dir, name))
str(data)

# Pathway gmt file obtained from
# https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
# https://yulab-smu.top/biomedical-knowledge-mining-book/useful-utilities.html?q=gene_set#read-gmt

pathway <- paste0(path_dir, "/c5.go.bp.v2023.1.Hs.entrez.gmt")

## Gene List set up  ------------------------------------------------------------
# Get enterez ids for pathway enrichment
genes <- data$hgnc_symbol

# Get gene info database
src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene", overwrite = TRUE)
# Map gene symbol to enterez ids
genes_ez <- Organism.dplyr::select(src, 
                                   keys = genes,
                                   columns = c("symbol","entrez"),
                                   keytype = "symbol")
genes_ez <- distinct(genes_ez)

# add back to orignal df
genes_ez <- rename(genes_ez, hgnc_symbol = symbol)

data <- left_join(data, genes_ez, by = "hgnc_symbol")

# Over representation analysis -------------------------------------------------
## Gene ontology enrichment -----------------------------------------------------

ego <- enrichGO(gene = genes_ez$entrez,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", readable = TRUE)
View(ego@result)

# Disease ontology
edo <- enrichDO(gene = genes_ez$entrez,
                ont = "DO",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", readable = TRUE)

View(edo@result)

## Visualisation ---------------------------------------------------------------

dotplot(ego, showCategory=15)
dotplot(edo, showCategory=15)

# Cluster mapping
p1 <- cnetplot(ego)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(ego, categorySize="pvalue")
p3 <- cnetplot(ego, circular = TRUE, colorEdge = TRUE) 

p1
p2
p3

## C7/C5 pahway enrichment ---------------------------------------------------------

msig <- read.gmt(pathway)
msig_ego <- enricher(genes_ez$entrez, TERM2GENE=msig)
msig_sig <- subset(msig_ego@result, p.adjust <0.05)

View(msig_ego@result)

dotplot(msig_ego, showCategory=15)
# Note C5 is the same db I believe used in the default function, but here you can customise which gmt database you use for enrichment

# TODO tidy plot labels

# Gene Set Enrichment Analysis -------------------------------------------------
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html#gsea-algorithm

### Step 1: Prepare gene list -----------------------------
# https://yulab-smu.top/biomedical-knowledge-mining-book/faq.html?q=geneList#genelist

# Remove duplicate genes
data <- data[!duplicated(data$entrez), ]

## feature 1: numeric vector: t-statsitic, LogFC, or P
geneList = data[, "t"]

## feature 2: named vector, entrez_ids from prev section
names(geneList) = as.character(data[,"entrez"])

## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)
head(geneList)

## Step 2: Prepare gene set ------------------------------------
# read gmt file of choice

msig <- read.gmt(pathway)

gsea <-
  GSEA(
  geneList,
  exponent = 1,
  minGSSize = 5,
  maxGSSize = 250,
  eps = 1e-10,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  TERM2GENE = msig,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea"
)

View(gsea)

# TODO compare clusters
# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-comparecluster.html
# Take top 3 sig HBG related clusters and compare/visualise: cnetplot()