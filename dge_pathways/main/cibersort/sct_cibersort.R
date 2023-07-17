### R set up and Functions -----------------------------------------------------

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

# load data 

load(file = paste0(results_dir, "/sct_eqtl.RData"))

dge <- dge[row.names(dge$counts) %in% filtered_genes$Name, ]
dge$genes <- filtered_genes

tpm <- as.data.frame(tpm(dge$counts, dge$genes$size))
head(tpm)

tpm$hgnc_symbol <- dge$genes$hgnc_symbol
tpm$dup <- stri_duplicated(tpm$hgnc_symbol)
table(tpm$dup)

tpm <- tidylog::filter(tpm, stri_duplicated(hgnc_symbol)==FALSE) 
tpm <- tidylog::filter(tpm, str_detect(hgnc_symbol, "-")==FALSE)
tpm <- select(tpm,-dup)

# reorder cols
tpm <- tpm[, c(806,1:805)]
row.names(tpm) <- tpm$hgnc_symbol

write.table(tpm, file = paste0(results_dir, "/sct_tpm_rnaseq.txt"), row.names = F, sep = "\t", quote = F)

## Gene file --------------

genes_c <- select(trans_sig, hgnc_symbol)
genes_d <- select(cis_sig, hgnc_symbol)

genes_c <- rbind(genes_d, genes_c)

write.table(genes_c, file = paste0(results_dir, "/trans_subset.txt"), row.names = F, col.names = F, sep = "\t", quote = F)


# Results -------------------------------------

# Cell fractions

# Cell expression indv

tcells <- t(fread(file =paste0(results_dir, "/ciber/CIBERSORTxHiRes_Job1_TcellsCD4_Window16.txt")))
tcells <- as.data.frame(row_to_names(tcells, row_number = 1))

# filter na and 1 genes
tcells <- tcells %>% select(where(not_all_na))

# Filter no variance genes
tcells <- apply(tcells, MARGIN=c(1,2), as.numeric)
tcells <- tcells[,colSums(tcells)!=nrow(tcells)]

library(caret)
nzv <- nearZeroVar(tcells)
tcells2 <- tcells[,  -nzv]
dim(tcells2)

# 154 genes with variance

# DEG cell-specific ------------------------------------------------------------

geno_pheno <- as.data.frame(cbind(pheno$sct, tcells2))

geno_pheno <- geno_pheno %>% rename(sct=V1)
geno_pheno$sct <- as.factor(geno_pheno$sct)

comp <- list(c("0", "1"))

ggplot(geno_pheno, aes(x=sct, y=ICOS, fill=sct)) +
  geom_boxplot() +
  labs(x="SCT", title = "CD4 T-cells") +
  scale_fill_manual(values = c("#FEA82F", "#A72020")) +
  theme_light() + theme(legend.position = "none") +
  stat_compare_means(comparisons = comp, paired = FALSE, method = "wilcox.test")
