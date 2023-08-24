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
library(ReactomePA)

# Set up:
source("~/scripts/r-packages.R")
source("~/scripts/functions.R")

## Directories ------------------------------------------------------

wd <- "~/Documents/whi_sca/rna/results" # Where main data is saved
results_dir <- "~/Documents/whi_sca/rna/results/dge_enrich" # "~/Documents/whi_sca/rna/results/wgcna"
plot_dir <- "~/Documents/whi_sca/rna/results/dge_enrich" # "~/Documents/whi_sca/rna/plots/wgcna/updates/just_sct/no_bc" # Where to save output plots
path_dir <- "~/Downloads/public_data" # Download from msigDB
setwd(wd)
load(file.path(results_dir, "pathway_fgsea.RData"))

## Input ------------------------------------------------------
# Pathway gmt file obtained from
# https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
# https://yulab-smu.top/biomedical-knowledge-mining-book/useful-utilities.html?q=gene_set#read-gmt

# Gene ontology
c5_all <- paste0(path_dir, "/c5/c5.all.v2023.1.Hs.entrez.gmt")
c5_bp <- paste0(path_dir, "/c5/c5.go.bp.v2023.1.Hs.entrez.gmt")
c5_cc <- paste0(path_dir, "/c5/c5.go.cc.v2023.1.Hs.entrez.gmt")
c5_mf <- paste0(path_dir, "/c5/c5.go.mf.v2023.1.Hs.entrez.gmt")
#  Pathways
c2_all <- paste0(path_dir, "/c2/c2.all.v2023.1.Hs.entrez.gmt")
#  Immune
c7_all <- paste0(path_dir, "/c7/c7.immunesigdb.v2023.1.Hs.entrez.gmt")


# Gene info database
src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene", overwrite = TRUE)

# OPTION WGCNA or DGE output
# Input data should be a list of gene symbols or ensembl/enterez ids to match against db
# name <- "gene_list_violet.csv" # WGCNA
# name <- "top_sig_dge_fc0.csv" # DGE filtered results
name <- "top_tab.csv" # DGE all results
data <- read.csv(file = file.path(wd, name)) # results_dirx
str(data)

## Gene List set up  ------------------------------------------------------------
# Get enterez ids for pathway enrichment
data$hgnc_symbol <- data$Description
genes <- data$hgnc_symbol

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

# Just perform with top significant genes:
top_data <- read.csv(file = file.path(wd, name)) 
top_data <- filter(top_data, adj.P.Val < 0.05)
top_data$hgnc_symbol <- top_data$Description
genes <- top_data$hgnc_symbol

# Map gene symbol to enterez ids
genes_ez2 <- Organism.dplyr::select(src, 
                                   keys = genes,
                                   columns = c("symbol","entrez"),
                                   keytype = "symbol")
genes_ez2 <- distinct(genes_ez2) # 136 genes

# add back to orignal df
genes_ez2 <- rename(genes_ez2, hgnc_symbol = symbol)
top_data <- left_join(top_data, genes_ez2, by = "hgnc_symbol")

## Gene ontology enrichment
ego <- enrichGO(gene = genes_ez2$entrez,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pvalueCutoff = 0.05,
                minGSSize = 3,
                maxGSSize = nrow(top_data),
                pAdjustMethod = "BH", readable = TRUE)

View(ego@result)
dotplot(ego)

# Disease ontology
edo <- enrichDO(gene = genes_ez2$entrez,
                ont = "DO",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", readable = TRUE)

View(edo@result)

## ORA Cluster Visualisation ---------------------------------------------------------------

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

## C7/C5 pathway enrichment ---------------------------------------------------------

path <- enrichPathway(gene=genes_ez2$entrez, pvalueCutoff = 0.05, readable=TRUE)
View(path@result)
dotplot(path)

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
geneList <- data[, "t"]
geneList <- data[, "logFC"]

## feature 2: named vector, entrez_ids from prev section
names(geneList) = as.character(data[,"hgnc_symbol"])
## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)
head(geneList)

## Step 2: Prepare gene set ----------------------------------------------------
# Read gmt file of choice

pathway <- c5_all
msig <- read.gmt(pathway)
# convert back into gene names
msig$entrez <- msig$gene
msig <- left_join(msig, genes_ez, by = "entrez")
msig <- msig %>% select(term, hgnc_symbol) %>% rename(gene = hgnc_symbol) %>%
  filter(!is.na(gene))

set.seed(1234)
gsea <-
  GSEA(
  geneList,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 0,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  TERM2GENE = msig,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea"
)

# Tidy results
View(gsea@result)
gsea@result$p.adjust <- signif(gsea@result$p.adjust, digits = 3)
gsea@result$Description <- str_remove(gsea@result$Description, ".*?_")

# Plot
dotplot(gsea, showCategory = 10) 

## C2 Curated gene sets e.g kegg -----------------------------------------------
c2 <- gsea
c2 <- pairwise_termsim(c2)
dotplot(c2, showCategory = 20)
treeplot(c2, showCategory = 20, nWords=3, nCluster = 5, hclust_method = "ward.D2")

# Filter results 
path_select <- c2$Description[1:25]
path_select <- path_select[!grepl("UV|INFLUENZA|CANCER", path_select, ignore.case = TRUE)]
path_select
# Filter similar terms
dotplot(react, showCategory = path_select)

# Reactome
react <- c2
react@result <- filter(c2@result, str_detect(c2@result$ID, "REACT"))
View(react@result)
dotplot(react, showCategory = 20)
# Filter results 
path_select <- react$Description[1:25]
path_select <- path_select[!grepl("SARS|INFLUENZA|ORC1", path_select, ignore.case = TRUE)]
path_select
# Filter similar terms
path_select <- path_select[!grepl("SARS|INFLUENZA|ORC1", path_select, ignore.case = TRUE)]

dotplot(react, showCategory = path_select)
treeplot(react, showCategory = path_select, nWords=3, nCluster = 4, hclust_method = "ward.D2")

# Kegg
kegg <- c2
kegg@result <- filter(c2@result, str_detect(c2@result$ID, "KEGG"))
dotplot(kegg, showCategory = 20) 
treeplot(kegg, showCategory = 15, nWords=3, nCluster = 10, hclust_method = "ward.D2")

# kidney related
kid <- c2
kid@result <- filter(c2@result, str_detect(c2@result$Description, "KIDNEY"))


## C5 ontology gene sets -------------------------------------------------------
c5 <- pairwise_termsim(c5)

# Subset
bp <- c5
bp@result <- filter(bp@result, str_detect(bp@result$ID, "BP"))
bp <- pairwise_termsim(bp)
dotplot(bp, showCategory = 10)
treeplot(bp, showCategory = 35, nWords=3, nCluster = 5, hclust_method = "ward.D2")

# Subset
path_select <- bp$Description[!grepl("HINDBRAIN|NEGATIVE|REGULATION", bp$Description, ignore.case = TRUE)]
path_select
treeplot(bp, showCategory = path_select, nWords=3, nCluster = 6, hclust_method = "ward.D2")
MF <- c5
MF@result <- filter(MF@result, str_detect(MF@result$ID, "MF"))
dotplot(MF, showCategory = 10)

cc <- c5
cc@result <- filter(cc@result, str_detect(cc@result$ID, "CC"))
dotplot(cc, showCategory = 10)

HP <- c5
HP@result <- filter(HP@result, str_detect(HP@result$ID, "HP"))
dotplot(HP, showCategory = 10)


## C7 immune -------------------------------------------------------------------

c7 <- gsea
dotplot(c7, showCategory = 10)


# Network plots ----------------------------------------------------------------

bp <- pairwise_termsim(bp)
MF <- pairwise_termsim(MF)
HP <-  pairwise_termsim(HP)
react <- pairwise_termsim(react)

treeplot(bp, showCategory = 20, nWords=3, nCluster = 5, hclust_method = "ward.D2")
treeplot(MF, showCategory = 20)
treeplot(HP, showCategory = 20)

# rm sars
react_filt <- react
react_filt@result <- filter(react_filt@result, !str_detect(react_filt@result$Description, "SARS"))

treeplot(react_filt, showCategory = 20, nWords=3)

c2 <- pairwise_termsim(c2)
treeplot(c2, showCategory = 30)

## KEGG -------------------------------------------------------------------
names(geneList) = as.character(data[,"entrez"])
## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)
head(geneList)

mkk2 <- gseMKEGG(geneList = geneList,
                 organism = 'hsa',
                 pvalueCutoff = 1)
head(mkk2)
View(mkk2@result)

# Reactome ---------------------------------------------------------------

react2 <- gsePathway(geneList, 
                pvalueCutoff = 0.2,
                pAdjustMethod = "BH", 
                verbose = FALSE)
head(react2)
View(react2@result)


## Save ------------------------------------------------------------------------
write.csv(c7, file = file.path(results_dir, "c7_dge_fgsea_t.csv"), row.names = F, quote = F)
write.csv(c5, file = file.path(results_dir, "c5_dge_fgsea_t.csv"), row.names = F, quote = F)
write.csv(c2, file = file.path(results_dir, "c2_dge_fgsea_t.csv"), row.names = F, quote = F)


save.image(file.path(results_dir, "pathway_fgsea.RData"))


# TODO compare clusters
# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-comparecluster.html
# Take top 3 sig HBG related clusters and compare/visualise: cnetplot()

