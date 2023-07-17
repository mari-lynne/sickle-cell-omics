# Plan
# Use top dge genes and see if they stratify with egfr/cell counts or can predict it

# Set up data ------------------------------------------------------------------

source("~/scripts/r-packages.R")
source("~/scripts/functions.R")

# Options and Directories ------------------------------------------------------

wd <- c("~/Documents/whi_sca/rna") # Where main data is saved
plot_dir <- c("~/Documents/whi_sca/rna/plots") # Where to save output plots
id_dir <- ("~/Documents/whi_sca/rna/ids")
meta_dir <- c("~/Documents/whi_sca/rna/meta")
results_dir <- c("~/Documents/whi_sca/rna/results")
setwd(wd)

# DGE data
load(file = paste0(results_dir, "/sct_dge.RData"))

# Kidney data
whi <-
  read.csv(file = paste0(meta_dir, "/WHI_HTNKidney_2023-02-27.csv")) %>% 
  rename(subject_id = subjectID)

# Subset sct data 
sct <- dge[, (dge$samples$group == "1")]  

# WHI pheno data
pheno_sct <- sct$samples
whi_sct <- filter(whi, subject_id %in% pheno_sct$subject_id)
whi_sct <- left_join(pheno_sct, whi_sct, by = "subject_id") %>% select(where(not_all_na))
# CKD
table(whi_sct$ckd)
# 791 nCKD, 46 = CKD (healthy) - 140:12 in SCT
# Too small sample size for disease classification, so focus on eGFR

# eGFR -------------------------------------------------------------------------


# egfr-ckd epi
# egfr-ckd modified MDRD

whi <- filter(whi, subject_id %in% pheno$subject_id)
whi <- left_join(pheno, whi, by = "subject_id") %>% select(where(not_all_na))

whi <- whi %>% rename(eGFR = eGFRCKDEpi)

p <- whi %>%
  ggplot(aes(x= eGFR, fill = sct)) +
  geom_density(alpha = 0.4) + 
  scale_fill_manual(values = c("#A72020", "#FEA82F")) +
  xlab("eGFR CKD-Epi")  +
  theme_light() +
  facet_wrap(~sct) 

p

png(paste(plot_dir,"egfr_ckdepi_sct.png", sep ="/"), width = 1880, height = 1680, res = 300)
p
dev.off()

# Slight binomial distribution
# could split into eGFR high versus low

whi %>% filter(sct == "1") %>%
  ggplot(aes(x= eGFR, fill = sct)) +
  geom_density(alpha = 0.4) + 
  scale_fill_manual(values = c("#FEA82F")) +
  xlab("eGFR CKD-Epi")  +
  theme_light()

# Notes
# Machine learning predict eGFR ------------------------------------------------

# Data we need
# eGFR phenotype
# Names of top DGE genes FDR < 0.01 (start with this for now)
# Top DGE genes overlap with SCD top DGEs, FDR < 0.01

# Then subset to only include top DGEs - raw counts??

# Recursive feature selection
# use log2cpm from v - try to intergrate weights at some point



# DGE data ---------------------------------------------------------------------

top <- toptab %>% filter(adj.P.Val < 0.05) # 113 genes

idx <- (v$genes$Description %in% top$gene_name)
idx2 <- (v$targets$group == 1)
sub <- v[idx, idx2]
dim(sub)

pheno <- sub$targets
pheno <- left_join(pheno, whi, by = "subject_id") %>% select(where(not_all_na))
pheno <- pheno %>% select(!ends_with(".x"))
colnames(pheno) <- str_remove_all(colnames(pheno), ".y")

genes <- t(sub$genes)
counts2 <- as.data.frame(t(sub$E))
colnames(counts2) <- genes[2,]

pheno_exprs <- cbind(pheno, counts2)

library(caret)

# Feature selection, selects genes which best correlate with eGFR
# Built in feature selection model - glmnet

# Data pre processing ----------------------------------------------------------

counts2 <- as.matrix(counts2)
# Remove zero variance

nzv <- nearZeroVar(counts2, saveMetrics= TRUE)
# no zero variance genes

# correlated predictors check

cor_genes <-  cor(counts2)
high_cor <- findCorrelation(cor_genes, cutoff = .75) # 296 samples
counts3 <- counts2[,-high_cor]
dim(counts3) # 145 genes

# Data splitting ---------------------------------------------------------------

# Make DF containing gene expression and egfr vars
pheno_exprs <- cbind(pheno, counts3)

set.seed(3456)
trainIndex <- createDataPartition(pheno_exprs$eGFR, 
                                  p = 0.75, 
                                  list = FALSE,
                                  times = 1)
head(trainIndex)

# get training set
datTrain <- pheno_exprs[trainIndex, ]
# get test set
datTest <- pheno_exprs[-trainIndex, ]

# Train data ---------------------------------------------------------------
mod <- caret::train(eGFR ~., data = datTrain, method = "glmnet", preProcess("center", "scale"))

mod

# linear regression -----------------------------------------------------------

mod <- lm(eGFR ~ counts3 + PC1 + PC2 + PC3 + PC4 + PC5 +PC6 +PC8 +PC9 + PC10,
          data = pheno_exprs)
summary(mod)

> design <- model.matrix( ~TG)
> y <- voom(dge, design)
> fit <- lmFit(y, design)
> fit <- eBayes(fit)
> diffTable2 <- topTable(fit, coef="TG", number=15)







# Just TAL1 --------------------------------------------------------------------

tal1 %>% filter(sct == "1") %>%
  ggplot(aes(x = sct, y = TAL1_counts)) +
  geom_boxplot()+
  theme_light() +
  scale_y_continuous(trans='log10') +
  labs(y = "TAL1 expression (counts)\n", x = "SCT status") +
  stat_compare_means(comparisons = comp)

