# wgcna

# Set up:

source("~/scripts/r-packages.R")
source("~/scripts/functions.R")
# source("~/scripts/color-palettes.R")

# Extra packages 

library(magrittr)
library(DESeq2)

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
# pheno <- filter(pheno, sct == "1" & flagged_seq_qc != "Yes") # 144 SCT participants
pheno <- filter(pheno, flagged_seq_qc != "Yes")

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
# summary_phen<- summarize_phenotypes(pheno = covar,
                                    # categorical_variables = c("smoking"),
                                    # numeric_variables = c("age","eGFR","eGFRMDRD"),
                                    # strata = "smoking")

### Filter SCT ----------------------------------------------------------------
## Genes and counts ------------------------------------------------------------

# Filter data to just include sct genotyped samples, and remove name/descript cols 
sub <- as.matrix(data[,colnames(data)%in% covar$rnaseq_ids])
# TODO (turn into an option)
sct_ids <- covar[covar$sct == 1, rnaseq_ids] # 144
sct_covar <- covar[covar$sct == 1, ]

sub <- sub[, colnames(sub) %in% sct_ids]

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
dge <- DGEList(counts=sub, samples=sct_covar, genes=filter_genes) # covar for all samps

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
keep.exprs <- filterByExpr(dge, group = dge$samples$sct)
dge <- dge[keep.exprs, ,keep.lib.sizes=FALSE]
dim(dge) # 13938  144

## Batch correction ------------------------------------------------------------

#TODO add PCA plots
exprs <- as.matrix(dge$counts)
pData <- dge$samples
pData$sct <- as.factor(pData$sct)
pData <- select(pData, plate, sct) # TODO test with and without sct covar

combat <- sva::ComBat_seq(counts=exprs, batch=pData$plate, full_mod=TRUE) # group=pData$sct
row.names(combat) <- dge$genes$Name

dim(combat) # rows = genes, columns = samples, needs to be transposed for wgcna format rows = samples, columns = genes

# WGCNA set up -----------------------------------------------------------------

### Normalise -------------------------------------------------------------------

combat_norm <- vst(combat)

# Update DGE object
dge_bc <- DGEList(counts=combat_norm, samples=dge$samples, genes=dge$genes)

## HBF Phenotype ---------------------------------------------------------------

# Subset data for just Hemoglobin expression data
idx <- (dge_bc$genes$hgnc_symbol %in% c("HBB", "HBG1", "HBG2")) 
sub <- dge_bc[idx,]
dim(sub)

# Make table with HBB and HBG counts, ratios and SCT status for plotting
pheno <- sub$samples
genes <- t(sub$genes)
counts2 <- as.data.frame(t(sub$counts))
colnames(counts2) <- genes[2,]

# Calculate ratios
# For now just aggregate HBG1/2
counts2$HBG <- counts2$HBG1 + counts2$HBG2
counts2$Ratio_HBG <- counts2$HBB / counts2$HBG
pheno_exprs <- cbind(pheno, counts2)
dim(pheno_exprs)

### WGCNA format ---------------------------------------------------------------
# wgcna format rows = samples, columns = genes

gene_matrix <- t(dge_bc$counts)

# Clinical phenotypes + cell counts --------------------------------------------
clin <- clean_names(read.csv(file = paste0(meta_dir, "/WHI_inflamation_2023-02-27.csv")))
whi_clin <- filter(clin, subject_id %in% pheno$subject_id)
# pheno_exprs
pheno_exprs <- left_join(pheno_exprs, whi_clin, by= "subject_id")
pheno_exprs <- pheno_exprs %>% select(-ends_with(".y"))
colnames(pheno_exprs) <- str_remove_all(colnames(pheno_exprs), ".x")
pheno_exprs <- pheno_exprs %>% select(-c("as311", "as315", "baa23"))
# write.csv(pheno_exprs, file = paste0(meta_dir, "/sct_only_study_clin_hbb_pheno.csv"))

# 0) exlude sampe outliters ----------------------------------------------------

sampleTree = hclust(dist(gene_matrix), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(10,7)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# TODO redo but need to add in the fact that some are SCT so cluster seprately and look for outliers separetly

# Plot a line to show the cut
abline(h = 70, col = "red");
# Determine cluster under the line # 80 # minSize 10
clust = cutreeStatic(sampleTree, cutHeight = 70, minSize = 10)
table(clust)
# clust 0 are outliers

keepSamples = (clust!= 0) # 19 outliers
gene_matrix = gene_matrix[keepSamples, ]
nGenes = ncol(gene_matrix)
nSamples = nrow(gene_matrix)
# datTraits = pheno_exprs[keepSamples, ]
# pick trait Ratio_HBG
datTraits = pheno_exprs[keepSamples, "Ratio_HBG"]

# Re-cluster samples
sampleTree2 = hclust(dist(gene_matrix), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")


# 1)  Pick soft threshold ----------------------------------------------------------

# Decides on what threshold a strong correlation will be defined as, later in the network
# So a higher power i.e 10 means that the scale free model tends to have 0.8 R2 which is hgh

library(WGCNA)
# allowWGCNAThreads()

powers <- c(1:14)
sft <- pickSoftThreshold(gene_matrix, powerVector = powers)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# function suggests 6 but I am going to use 8 based on inspection of scale free top model fit and mean connectivity

### Portal-1 -------------------------------------------------------------
rm(covar, filter_genes)
#save.image(file = paste0(results_dir, "/1_sct_only_wgcna.RData"))
# load(file = paste0(results_dir, "/1_sct_only_wgcna.RData)
#TODO open_portal function

# 2) Perform signed network analysis -------------------------------
allowWGCNAThreads()
# Unsupervised clustering analysis
# Create a signed network using the chosen soft-thresholding power
signed_net <- blockwiseModules(
  gene_matrix,
  power = 9, # sft$powerEstimate, # Use the estimated power value from previous step 
  weights = NULL,
  TOMType = "signed",
  networkType = "unsigned",
  corType = "pearson",
  saveTOMs = FALSE,
  verbose = 3
)

rm(combat, combat_norm, counts2, dge, ensembl)

### Portal-2----------------------------------------------
# save.image(file = "wgcna.RData")
# save.image(file = paste0(results_dir, "/2_sct_only_wgcna.RData"))
# load(file = paste0(results_dir, "/2_sct_only_wgcna.RData")

load(file = paste0(results_dir, "/wgcna.RData"))

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(signed_net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(signed_net$dendrograms[[1]], mergedColors[signed_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


moduleLabels = signed_net$colors
moduleColors = labels2colors(signed_net$colors)
MEs = signed_net$MEs;
geneTree = signed_net$dendrograms[[1]];

save(MEs, moduleLabels, moduleColors, geneTree,
     file = "sct-all-networkConstruction.RData")

# 3 Identify modules (clusters) of co-expressed genes: -------------------------

## Select phenotypes

datTraits = pheno_exprs[keepSamples, ]
datTraits <- datTraits %>% select(eGFR, bmi_t0, Ratio_HBG, hgb, rbc_count, rbc_dist_width, platelet_count, wbc, lymphocytes)

datTraits <- datTraits %>% select(Ratio_HBG, lymphocytes)

# Define numbers of genes and samples
nGenes = ncol(gene_matrix);
nSamples = nrow(gene_matrix);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(gene_matrix, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
### Module correlation Plot ------------------------------
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
ylabs <- str_remove_all(names(MEs), "ME")
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = ylabs,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = FALSE,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


library(ggplot2)
# Calculate the number of genes per module
moduleGenes <- table(moduleColors)

# Merge module colors with the number of genes
mergedData <- data.frame(ModuleColor = names(moduleGenes), NumGenes = as.numeric(moduleGenes))

# Initialize the plot and specify aesthetics
plot <- ggplot(mergedData, aes(x = ModuleColor, y = NumGenes))

# Add bar plot layer
plot <- plot + geom_bar(stat = "identity", fill = "blue") + labs(x = "Module Color", y = "Number of Genes", title = "Number of Genes per Module")

# Display the plot
plot




## Focus on HBF ------------------------------------------

# Define variable weight containing the weight column of datTrait
HBF = as.data.frame(datTraits$Ratio_HBG);
names(HBF) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(gene_matrix, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(gene_matrix, HBF, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(HBF), sep="");
names(GSPvalue) = paste("p.GS.", names(HBF), sep="");


# Investigate genes in modules
names(gene_matrix)
mod_col = c("sienna3")
gene_interest <- colnames(gene_matrix)[moduleColors== mod_col]
# Filter .11
gene_interest <- sub("\\.[0-9]+$", "", gene_interest)

# Get gene names

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(
  attributes=c('ensembl_gene_id', 'hgnc_symbol','chromosome_name','start_position','end_position'),
  mart = ensembl)

gene_interest <- filter(genes, ensembl_gene_id %in% gene_interest)

write.csv(gene_interest, file = paste0(results_dir, "/gene_interest_sienna3.csv"))

# compare against toptab
gene_interest <- as.data.frame(gene_interest)

toptab <- read.csv(file = "/home/mari/Documents/whi_sca/rna/results/toptab.csv")
toptab <- filter(toptab, (adj.P.Val < 0.05 & (logFC < -0.1 | logFC > 0.1 ))) %>% select(Description)

# 47 genes from trait are in this module
overlap <- gene_interest[gene_interest$hgnc_symbol %in% toptab$Description, ]

dim(overlap)







# Extra ------------

# 3 Identify modules (clusters) of co-expressed genes: -------------------------

# Module assignment

# Choose a minimum module size (number of genes)
min_module_size <- 30

# Merge similar modules and define the module colors
merged_modules <- mergeCloseModules(signed_net, cutHeight = 0.25, verbose = 3, colors = module_colors)
module_colors <- labels2colors(signed_net$colors)

# Assign modules to genes
module_assignment <- merged_modules$colors
names(module_assignment) <- colnames(gene_matrix)

# 4 Identify modules significantly associated with HBF expression: -------------

phenotype_data <- pheno

# Convert HBF expression vector to numeric
hbf_expression <- as.numeric(phenotype_data$HBF_Expression)

# Perform module-trait relationship analysis
module_trait_cor <- cor(hbf_expression, module_assignment, use = "p")
module_trait_pvalue <- corPvalueStudent(module_trait_cor, nSamples = length(hbf_expression))

# Select modules with significant correlation and adjust p-values for multiple testing
module_trait_threshold <- 0.05  # Adjust as needed
sig_modules <- module_trait_pvalue < module_trait_threshold

# 5)  Visualise the analysis -------------------------------------------

# Plot the network heatmap
plotHeatmap(
  signed_net,
  main = "WGCNA Heatmap",
  moduleColors = module_colors
)

# Plot the module-trait relationship
plotModuleTraitHeatmap(
  gene_matrix,
  moduleColors = sig_module_colors,
  traitData = hbf_expression,
  main = "Module-Trait Relationship Heatmap"
)

sig_module_colors <- module_colors[sig_modules]

