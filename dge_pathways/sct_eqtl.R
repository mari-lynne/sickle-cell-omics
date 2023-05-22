## 1) Pre-processing

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

# Data Set up ------------------------------------------------------------------

# RNA seq
data <- read_omic(name = data_name, wd = wd)

## Phenotype and meta data -----------------------------------------------------
# (see sct_id_pheno_match.R for more)

pheno <- fread(paste(meta_dir, pheno_name, sep = "/"))
pheno <- filter(pheno, flagged_seq_qc != "Yes") # 805

# Filter data to just include sct genotyped samples 
counts <- data[,colnames(data)%in% pheno$rnaseq_ids]

# Update counts
counts <- as.matrix(counts)
row.names(counts) <- gene_info$Name

# Sum of mapped reads per sample
lsize <- colSums(counts)
pheno$lib.size_og <- lsize # save in pheno 

### Add PCA data ---------------------------------------------------------------

# See pca.R for PCA calculations
pc_top <- fread(file = paste0(results_dir, "/pc_top.txt"))
pc_top <- pc_top %>% rename(rnaseq_ids = V1)

covar <- left_join(pheno, pc_top)

# Update other covars, clean names, check factors
covar <- covar %>%
  mutate(plate = coalesce(covar$pick_plate1, covar$pick_plate_2_repeat)) %>%
  select(-c("pick_plate1", "pick_plate_2_repeat"))

covar$plate <- as.factor(covar$plate)
covar <- covar %>% mutate(sct = ifelse(sct == 1, "SCT", "Control"))
covar$sct <- as.factor(covar$sct)

# Recode smoking never missing vs current and past
covar <- covar %>%
  mutate(smoking = ifelse
         (smoking == "Missing", "Never", covar$smoking))

### WHI clinical covars --------------------------------------------------------

# Add whi data
whi <-
  read.csv(file = paste0(meta_dir, "/WHI_HTNKidney_2023-02-27.csv")) %>% 
  rename(subject_id = subjectID)

whi_sct <- filter(whi, subject_id %in% pheno$subject_id)
covar <- left_join(covar, whi_sct, by= "subject_id") %>% select(where(not_all_na))

covar <- rename(covar, eGFR = eGFRCKDEpi)
covar <- covar %>% select("rnaseq_ids","commonid", "subject_id","bmi_t0", "smoking","age", "lib.size_og", "sct",
                           "PC1", "PC2","PC3","PC4", "PC5","PC6","PC7","PC8", "PC9", "PC10", "PC11", "PC12","PC13", "PC14",
                           "plate", "screat","sbp","ace","diet_poteassium", "diet_sodium",
                           "eGFR", "eGFRMDRD","ckd", "ageckd","AnyHTNmeds")

## Add gene info ---------------------------------------------------------------

gene_info <- data[, 1:2]
gene_info <- rename(gene_info, hgnc_symbol = Description)

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

genes <- getBM(
  attributes=c('ensembl_gene_id', 'hgnc_symbol','chromosome_name','start_position','end_position'),
  mart = ensembl)

genes$size = genes$end_position - genes$start_position

# Filter genes in counts table with missing gene info

filter_genes <- 
  left_join(gene_info %>% group_by(hgnc_symbol) %>% mutate(id = row_number()),
            genes %>% group_by(hgnc_symbol) %>% mutate(id = row_number()), 
            by = c("hgnc_symbol", "id")) %>%
            filter(!is.na(ensembl_gene_id)) %>% select(-id)

counts <- counts[row.names(counts) %in% filter_genes$Name, ]

### Form DGE list --------------------------------------------------------------

dge <- DGEList(counts=counts, samples=covar, genes=filter_genes, group = covar$sct)

# QC checks --------------------------------------------------------------------
lcpm <- cpm(dge$counts, log=TRUE)

L <- mean(dge$samples$lib.size) * 1e-6
M <- median(dge$samples$lib.size) * 1e-6
c(L, M) # 56.6 CPM - median

# Filter zero counts across all samples
zero_counts <- rowSums(dge$counts==0)==ncol(dge$counts)
dge <- dge[!zero_counts, ,keep.lib.sizes=FALSE]
dim(dge$counts)

# Filter low count
keep.exprs <- filterByExpr(dge)
dge <- dge[keep.exprs, ,keep.lib.sizes=FALSE]
dim(dge) # 15970   805


### MDS plots ------------------------------------------------------------------

# Filtered count data
lcpm2 <- cpm(dge, log = TRUE)
# Set up plate as group/factor to plot
group <- as.factor(dge$samples$plate) # Parse to labels
# Make colour vector using group levels
col.group <- group
# Viridis palette function
colours <- viridis_pal(alpha = 1, begin = 0, end = 1, direction = 1, option = "C")
levels(col.group) <- colours(nlevels(col.group))
# Recode levels into character vector of colors matching length of data
pal <- as.character(col.group) 
# Plot MDS, use PCA "common" option is faster
plotMDS(lcpm2, labels = group, col = pal, gene.selection = "common")

# Check principle components aren't capturing plate variation later
plotMDS(lcpm2, labels = group, col = pal, dim.plot = c(3,2), gene.selection = "common")

# SCT PCA
group <- as.factor(dge$samples$sct)
col.group <- group
colours <- viridis_pal(alpha = 1, begin = 0.6, end = 0, direction = 1, option = "A")
levels(col.group) <- colours(nlevels(col.group))
pal <- as.character(col.group) 
# Plot MDS, use PCA "common" option is faster
plotMDS(lcpm2, labels = group, col = pal, gene.selection = "common")
plotMDS(lcpm2, labels = group, col = pal, dim.plot = c(5,6), gene.selection = "common")


## Normalisation----------------------------------------------------------------

dge <- calcNormFactors(dge, method = "TMM")

# TODO fix box plot code

# Design Matrix ---------------------------------------------------------------
# Make var of interest a factor

design <-
  model.matrix(~0+ sct + age + plate + bmi_t0 + smoking +
                 PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10 +PC11 +PC12 +PC13 +PC14,
               data = dge$samples)

View(design)

# Limma-voom -------------------------------------------------------------------

pdf(paste(plot_dir,"filter_voom.pdf",sep ="/"))
v <- voom(dge, design, plot = TRUE)
dev.off()

pdf(paste(plot_dir,"filter_normcounts.pdf",sep ="/"))
boxplot(norm.expr, 
        col=col.group,
        main="Distribution of normalised counts",
        xlab="",
        ylab="log2 normalised expression",
        las=2,cex.axis=0.8)
dev.off()

# Fit a regression model to counts, with stabalised variance. Parameters supplied by design
vfit <- lmFit(v, design)

# Contrast matrix --------------------------------------------------------------

# Define contrasts to make
contr.mat <- makeContrasts(sct = sctSCT - sctControl, levels = colnames(design))
vfit <- contrasts.fit(vfit, contrasts=contr.mat)
efit <- eBayes(vfit)

pdf(paste(plot_dir,"filter_meanvar.pdf",sep ="/"))
plotSA(efit)
dev.off()

sct_top <- topTable(efit, coef = 1, number = Inf) # Adj = BH method with 0.05/5% FDR

trans_sig <- sct_top %>% filter(adj.P.Val <0.05) %>% select(-Name)

write.csv(trans_sig, file  = paste0(results_dir, "/transEQTL.csv"), row.names = FALSE)


# Volcanoes ---------------------------------------------------------------------

toptab <- sct_top
toptab$gene_name <- toptab$hgnc_symbol

plot_vol(toptab, lab = 'top', title = "DEG Sickle Cell Trait (rs334)",
         FC = 0.25, alpha = 0.6,
         colours = c("orange", "#ef5a3a", "#a50000", "#800000"))

plot_vol()

# Cis-eQTL genes ---------------------------------------------------------------

# Filter genes in sct_top for those in 1 MB of rs334
# Chromosome	11, Position	5227002, Gene	HBB]
start = 5227002-1e6
end = 5227002+1e6

cis <- sct_top %>%
  filter(chromosome_name == "11") %>%
  filter(start_position %in% seq(start,5227002)|end_position %in% seq(5227002, end))

cis <-
  cis %>% select(-adj.P.Val) %>%
  mutate(adj.P.Val = p.adjust(P.Value, method = "BH")) %>% select(-Name)

cis_sig <- filter(cis, adj.P.Val < 0.05)

write.csv(cis, file  = paste0(results_dir, "/cisEQTL.csv"), row.names = FALSE)

# https://maayanlab.cloud/Harmonizome/gene_set/sickle+cell+anemia/GWASdb+SNP-Disease+Associations


# Proteomics  ------------------------------------------------------------------

proteins <- fread(paste0(wd, "/meta/sct_proteins.txt"))
proteins$Protein <- str_trim(proteins$Protein)

# 23
proteins <- filter(sct_top, hgnc_symbol %in% proteins$Protein)

proteins <-
  proteins %>%
  mutate(adj.P.Val_2 = p.adjust(P.Value, method = "BH")) %>% select(-Name)

sig_proteins <- filter(proteins, P.Value < 0.05)

sig_proteins$Dataset <- rep("protein", nrow(sig_proteins))

# Meth -------------------------------------------------------------------------

meth <- read.csv(file = "sct_meth_genes.csv")
meth$Gene <- str_trim(meth$Gene)

meth <- filter(sct_top, hgnc_symbol %in% meth$Gene)

meth <-
  meth %>%
  mutate(adj.P.Val_2 = p.adjust(P.Value, method = "BH")) %>% select(-Name)

sig_meth <- filter(meth, P.Value < 0.05)

sig_meth$Dataset <- rep("meth", nrow(sig_meth))

meth_pro <- bind_rows(sig_meth, sig_proteins)

write.csv(meth_pro, file = paste0(results_dir, "/protein_methylation_overlap.csv"), row.names = FALSE)

### Save data -----
save.image(file = paste0(results_dir, "/sct_eqtl.RData"))


### JHS overlap ---------------------

jhs_cis <- fread("rs334_cis_JHS.txt")
jhs_cis <- filter(jhs_cis, `p-value` < 0.05) # 2 sig cis genes
jhs_whi_cis <- intersect(jhs_cis$gene, cis_sig$ensembl_gene_id)
jhs_whi_cis <- cis[cis$ensembl_gene_id %in% jhs_whi_cis, ]
# 3 have adj P < 0.05, 5 have P.val < 0.05
# filtered using non-adj p sig <0.05 JHS (x2 genes) and adj sig WHI (x5 genes) -> just UBQLN1 replciates

jhs_trans <- fread("rs334_trans_JHS.txt")
jhs_trans <- filter(jhs_trans, `p-value` < 0.05) # 856
jhs_whi_trans <- intersect(jhs_trans$gene, trans_sig$ensembl_gene_id)
jhs_whi_trans <- trans_sig[trans_sig$ensembl_gene_id %in% jhs_whi_trans, ]

# filtered using non-adj p sig <0.05 JHS (837 genes) and adj sig WHI (x299 genes) -> 18 genes replicate

trans2 <- sct_top %>% filter(P.Value <0.05) %>% select(-Name)
jhs_whi_trans2 <- intersect(jhs_trans$gene, trans2$ensembl_gene_id)
jhs_whi_trans2 <- trans2[trans2$ensembl_gene_id %in% jhs_whi_trans2, ]
# or with just filtering my results by non-adj P.val <0.05 = 114 replicate genes

### MESA overlap ----------------

mesa <- fread("MESA_qtl.csv")
mesa <- filter(mesa, pvalue <0.05) # 4335
mesa_whi_trans2 <- intersect(mesa$Gene_stable_ID, trans2$ensembl_gene_id)
mesa_whi_trans2 <- trans2[trans2$ensembl_gene_id %in% mesa_whi_trans2, ]
# 363 relaxed filtering

mesa_whi_trans <- intersect(mesa$Gene_stable_ID, trans_sig$ensembl_gene_id)
mesa_whi_trans <- trans_sig[trans_sig$ensembl_gene_id %in% mesa_whi_trans, ]
# 37 adj p filtering

# All overlap (trans)
reps <- intersect(jhs_whi_trans, mesa_whi_trans)

write.table(reps, file = paste0(results_dir, "/all_trans_qtl.csv"), row.names = F)
write.table(jhs_whi_trans2, file = paste0(results_dir, "/jhs_trans_qtl.csv"), row.names = F)
write.table(mesa_whi_trans2, file = paste0(results_dir, "/mesa_trans_qtl.csv"), row.names = F)

