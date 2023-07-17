# RNAseq for modelling linear traits

# Pipeline based on https://github.com/nkurniansyah/Olivia

# Set up:

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

# Load data --------------------------------------------------------------------  

# RNA seq
data <- read_omic(name = data_name, wd = wd)

## Pheno ------------------------------------------------------------------------

# Phenotype and meta data (see sct_id_pheno_match.R for more)
pheno <- fread(paste(meta_dir, pheno_name, sep = "/"))
pheno <- filter(pheno, sct == "1" & flagged_seq_qc != "Yes") # 144 SCT participants
# pheno <- filter(pheno, flagged_seq_qc != "Yes")

## PCA-covar -------------------------------------------------------------------------

# See pca.R
pc_top <- fread(file = paste0(results_dir, "/pc_top.txt"))
pc_top <- rename(pc_top, rnaseq_ids = V1)

covar <- left_join(pheno, pc_top, by = "rnaseq_ids")

covar <- covar %>%
  mutate(plate = coalesce(covar$pick_plate1, covar$pick_plate_2_repeat)) %>%
  select(-c("pick_plate1", "pick_plate_2_repeat"))

# Recode smoking never missing vs current and past
covar <- covar %>%
  mutate(smoking = ifelse
         (smoking == "Missing", "Never", covar$smoking))

# Update factors
covar$plate <- as.factor(covar$plate)

# Add whi data
whi <-
  read.csv(file = paste0(meta_dir, "/WHI_HTNKidney_2023-02-27.csv")) %>% 
  rename(subject_id = subjectID)

whi_sct <- filter(whi, subject_id %in% pheno$subject_id)
covar <- left_join(covar, whi_sct, by= "subject_id") %>% select(where(not_all_na))

covar <- rename(covar, eGFR = eGFRCKDEpi)

# Summarise pheno
summary_phen<- summarize_phenotypes(pheno = covar,
                                    categorical_variables = c("smoking"),
                                    numeric_variables = c("age","eGFR","eGFRMDRD"),
                                    strata = "smoking")
summary_phen

## Genes and counts ------------------------------------------------------------

# Filter data to just include sct genotyped samples, and remove name/descript cols 
sub <- as.matrix(data[,colnames(data)%in% covar$rnaseq_ids])

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
dge <- DGEList(counts=sub, samples=covar, genes=filter_genes)

# Filter genes missing ensembl and transcript length info 
dge <- dge[!is.na(dge$genes$ensembl_gene_id), ]

dim(dge)

# QC --------------------------------------------------------------------------

## Filter low counts

# Filter zero counts across all samples
zero_counts <- rowSums(dge$counts==0)==ncol(dge$counts)
dge <- dge[!zero_counts, ,keep.lib.sizes=FALSE]
dim(dge$counts)

# Filter low count
keep.exprs <- filterByExpr(dge)
dge <- dge[keep.exprs, ,keep.lib.sizes=FALSE]
dim(dge) # 13938  144

## Batch correction ------------------------------------------------------------

#TODO add PCA plots
exprs <- as.matrix(dge$counts)
pData <- dge$samples
pData <- select(pData, plate) # TODO test with and without sct covar

combat <- sva::ComBat_seq(counts=exprs, batch=pData$plate, full_mod=TRUE)
row.names(combat) <- dge$genes$Name

# Update DGE object
dge_bc <- DGEList(counts=combat, samples=dge$samples, genes=dge$genes)


## TPM -------------------------------------------------------------------------
options(scipen = 999)
tpm_bc <- as.data.frame(tpm(dge_bc$counts, dge_bc$genes$size))
row.names(tpm_bc) <- dge_bc$genes$Name

# Log2TPM
tpm_log <- apply(tpm_bc, 2, log2)
row.names(tpm_log) <- dge_bc$genes$Name

# Filter genes with low variance in respect to y
nzv <- nearZeroVar(tpm)
filteredDescr <- tpm[,  -nzv]
dim(filteredDescr)

# Modelling --------------------------------------------------------------------

## Data split -----------------------------------------------------------------

# Remove genes that do not correlate with eGFR
# Dependencies JDK (linux), rJAVA, RWeka, FSelector install_github("vqv/ggbiplot")
library(feseR)

# reorder so y var is last
# structure so n rows = samples, p cols = features
tpms <- as.data.frame(t(rbind(tpm_bc, dge_bc$samples$eGFR)))
y <- tpms[, ncol(tpms)]
y
# Split into test and training
library(caret)
set.seed(333)
trainIndex <- createDataPartition(tpms$`13939`, p = .8, 
                                  list = FALSE, 
                                  times = 1)

Train <- tpms[ trainIndex,]
Test  <- tpms[-trainIndex,]

test <- train(`13939`~.,data =Train, method ='glmStepAIC')

# Feature selection -----




## BOSO Modelling ---------------------------------------------------------------

tpm <- as.data.frame(t(rbind(dge_bc$samples$eGFR, tpm_log)))

# Split into test and training
library(caret)
set.seed(333)
trainIndex <- createDataPartition(tpm$V1, p = .8, 
                                  list = FALSE, 
                                  times = 1)

train <- tpm[ trainIndex,]
ytrain <- train$V1
train <- as.matrix(train)

test  <- tpm[-trainIndex,]
ytest <- test$V1
test <- as.matrix(test)

boso <- BOSO(x=train, y=ytrain,
             xval=test, yval=ytest,
             intercept = TRUE,
             seed=333)

data("sim.xy", package = "BOSO")

Xtr <- sim.xy[['high-5']]$x
Ytr <- sim.xy[['high-5']]$y
Xval <- sim.xy[['high-5']]$xval
Yval <- sim.xy[['high-5']]$yval

time <- Sys.time()
obj <- BOSO(x = Xtr, y = Ytr,
            xval = Xval, yval = Yval,
            IC = 'eBIC',
            nlambda=1000,
            intercept= 0,
            standardize = 0,
            Threads=4, timeLimit = 120, verbose = 3, 
            seed = 2021)
time <- as.numeric(Sys.time() - time)
