# Weighted Gene Correlation Network Analysis

# Set up:
source("~/scripts/r-packages.R")
source("~/scripts/functions.R")

# Options and Directories ------------------------------------------------------

wd <- c("~/Documents/whi_sca/rna") # Where main data is saved
plot_dir <- c("~/Documents/whi_sca/rna/plots") # Where to save output plots
meta_dir <- c("~/Documents/whi_sca/rna/meta")
results_dir <- c("~/Documents/whi_sca/rna/results")

setwd(wd)

# Load data --------------------------------------------------------------------  

# RNA seq
dge_bc <- readRDS(file = "dge_bc_all.rds")
# combat adjusted for plate, and also vst normalised after
covar <- read.csv(file = paste0(meta_dir, "/sct_all_covars_Jul23.csv"))

# Exclude hispanics for now
dge_bc <- dge_bc[, dge_bc$samples$ethnic==3]


### WGCNA format ---------------------------------------------------------------
# wgcna format rows = samples, columns = genes
# 
gene_matrix <- t(dge_bc$counts)

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

