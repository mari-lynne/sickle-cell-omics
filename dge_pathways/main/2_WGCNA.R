# Weighted Gene Correlation Network Analysis

# Set up:
source("~/scripts/r-packages.R")
source("~/scripts/functions.R")

# Options and Directories ------------------------------------------------------

wd <- c("~/Documents/whi_sca/rna") # Where main data is saved
plot_dir <- c("~/Documents/whi_sca/rna/plots/wgcna/updates/just_sct/no_bc") # Where to save output plots
meta_dir <- c("~/Documents/whi_sca/rna/meta")
results_dir <- c("~/Documents/whi_sca/rna/results")

setwd(wd)

# Load data --------------------------------------------------------------------  
# OPTION with our without bc
# RNA seq
dge <- readRDS(file = "dge_bc_black.rds")
dge <- readRDS(file = "dge_black.rds")

# combat adjusted for plate, and also vst normalised after
covar <- read.csv(file = paste0(meta_dir, "/sct_black_covars_hbg_Jul23.csv"))

# OPTION: Just SCT
dge <- dge[, dge$samples$sct==1]
# Option: remove manual outliers then redo sample clustering
outliers <- c("X941368")
dge <- dge[, dge$samples$rnaseq_ids %!in% outliers]

### WGCNA format ---------------------------------------------------------------
# wgcna format rows = samples, columns = genes
gene_matrix <- t(dge$counts)

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

# Cut first outlier in prev section, then redo # X941368

# TODO redo but need to add in the fact that some are SCT so cluster seprately and look for outliers separetly

# Plot a line to show the cut
abline(h = 80, col = "red");
# Determine cluster under the line # 80 # minSize 10
clust = cutreeStatic(sampleTree, cutHeight = 80, minSize = 10)
table(clust)
# clust 0 are outliers
keepSamples = (clust!= 0) # 5 outliers
gene_matrix = gene_matrix[keepSamples, ]
nGenes = ncol(gene_matrix)
nSamples = nrow(gene_matrix)

# OPTION: Pick trait 
datTraits = dge$samples[keepSamples, "Ratio_HBG"]

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
allowWGCNAThreads()

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(gene_matrix, powerVector = powers, verbose = 5)

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

print(sft$powerEstimate)
# function suggests 6 but I am going to use 8 based on inspection of scale free top model fit and mean connectivity
# 
#  Jul 14 redo of just sct, suggests 17?

### Portal-1 -------------------------------------------------------------

#save.image(file = paste0(results_dir, "/1_sct_only_wgcna.RData"))
# load(file = paste0(results_dir, "/1_sct_only_wgcna.RData)
#TODO open_portal function

#OPTION: custom sft power threshold
sft$powerEstimate <- 12

# 2) Perform signed network analysis -------------------------------
allowWGCNAThreads()
# Unsupervised clustering analysis
# Create a signed network using the chosen soft-thresholding power
signed_net <- blockwiseModules(
  gene_matrix,
  power = sft$powerEstimate, # Use the estimated power value from previous step
  weights = NULL,
  TOMType = "signed",
  networkType = "unsigned",
  corType = "pearson",
  saveTOMs = FALSE,
  verbose = 3
)

### Portal-2----------------------------------------------
# save.image(file = "wgcna.RData")
# load(file = paste0(results_dir, "/wgcna.RData"))
# save.image(file = paste0(results_dir, "/2_sct_only_wgcna.RData"))
load(file = paste0(results_dir, "/2_sct_only_wgcna.RData"))


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

## Portal 3 -----------------------------
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "3_sct-networkConstruction.RData")

# 3 Identify modules (clusters) of co-expressed genes: -------------------------

## Select phenotypes

datTraits <- dge$samples[keepSamples, ]
datTraits <- datTraits %>% select(Ratio_HBG, rbc_dist_width, neutrophil, tnf_alpha, eGFR, screat) # platelet_count bmi_t0  lymphocytes,

# datTraits <- datTraits %>% select(Ratio_HBG, wbc)

# Define numbers of genes and samples
nGenes = ncol(gene_matrix);
nSamples = nrow(gene_matrix);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(gene_matrix, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# Tidy
row.names(moduleTraitPvalue) <- str_remove_all(row.names(moduleTraitPvalue),"ME");
### Module correlation Plot ------------------------------

# Main plot
sizeGrWindow(10,6)
# OPTION: Pval labelling
# Will display correlations and their p-values
textMatrix <- ifelse(moduleTraitPvalue < 0.05, signif(moduleTraitPvalue, 1), "")
dim(textMatrix) <- dim(moduleTraitCor)

par(mar <- c(6, 8.5, 3, 3));
# OPTION ylab colour names or highlight
ylabs <- str_remove_all(names(MEs), "ME")
# Create the y-axis labels with "*" if the trait is significant
# merge data
# mutate
ylabs <- ifelse(moduleTraitPvalue[,"Ratio_HBG"] < 0.05, "*", "")
print(ifelse(moduleTraitPvalue[,"Ratio_HBG"] < 0.05, "sig", "NS"))

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = ylabs,
               colorLabels = FALSE,
               colors = hcl.colors(50, palette = "RdBu"),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("SCT Module-trait relationships"))

# violet,steelblue, orange
# steelblue, black, yellow, light cyan

# Calculate the number of genes per module
moduleGenes <- table(moduleColors)

# Merge module colors with the number of genes
mergedData <- data.frame(ModuleColor = names(moduleGenes), NumGenes = as.numeric(moduleGenes))

# Add bar plot layer
ngenes <-
  ggplot(mergedData, aes(x = ModuleColor, y = NumGenes, fill = ModuleColor)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = mergedData$ModuleColor) +
  labs(x = "Module Color", y = "Number of Genes",
       title = "Number of Genes per Module") +
         theme(legend.position = "none", axis.text.x = element_blank()) 


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

# Investigate genes in modules --------------------------

# Get a list of significant modules with trait of interest:
p_values <- moduleTraitPvalue[, "Ratio_HBG"]
# Extract the color names where p-value < 0.05
significant_colors <- rownames(moduleTraitPvalue)[p_values < 0.05]
# Print the significant color names
print(significant_colors)

# Get ensembl Mart and toptables
ensembl_mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
top_tab <- "/home/mari/Documents/whi_sca/rna/results/toptab.csv"

# Process clusters and overlapped genes
process_color_module(module_colors=significant_colors, gene_matrix=gene_matrix, ensembl_mart=ensembl_mart, results_dir=results_dir, top_table=top_tab)

# Portal-4: ----------------------------------------------------
# save.image(file = paste0(results_dir, "/4_wgcna_sct.RData"))
load(file = paste0(results_dir, "/4_wgcna_sct.RData"))


# FGSEA prep -------------------------------------------------------------------

# Need to list of gene names/ids plus the fold change and other stats from toptab
toptab <- read.csv(file = "/home/mari/Documents/whi_sca/rna/results/toptab.csv")
toptab <- rename(toptab, hgnc_symbol = gene_name)
geneList <- left_join(gene_interest, toptab, by = "hgnc_symbol")
geneList <- dplyr::select(geneList, c("hgnc_symbol", "ensembl_gene_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"))
write.csv(geneList, file = paste0(results_dir, "/gene_list_violet.csv"))

# Highlight plots ----------------------------------------

# Main plot
sizeGrWindow(10,6)
# OPTION: Pval labelling
# Will display correlations and their p-values
textMatrix <- ifelse(moduleTraitPvalue < 0.05, signif(moduleTraitPvalue, 1), "")
dim(textMatrix) <- dim(moduleTraitCor)

par(mfrow = c(1, 2));
# OPTION ylab colour names or highlight
ylabs <- str_remove_all(names(MEs), "ME")
# Create the y-axis labels with "*" if the trait is significant
# merge data
# mutate
ylabs <- ifelse(moduleTraitPvalue[,"Ratio_HBG"] < 0.05, "*", "")
print(ifelse(moduleTraitPvalue[,"Ratio_HBG"] < 0.05, "sig", "NS"))

png(paste0(plot_dir, "/module_trait_cor.png"), width = 7, height = 6, units = "in", res = 300)
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = ylabs,
               colorLabels = FALSE,
               colors = hcl.colors(50, palette = "RdBu"),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("SCT Module-trait relationships"))
dev.off()

# Highlight SCT cluster
png(paste0(plot_dir, "/ngenes_module.png"), width = 7, height = 6, units = "in", res = 300)
ngenes <- ngenes +
  annotate("text", x = "violet", y = moduleGenes["violet"], label = "*", color = "black", size = 8, vjust = 0)  +
  annotate("text", x = "steelblue", y = moduleGenes["steelblue"], label = "*", color = "black", size = 8, vjust = 0)  +
  annotate("text", x = "orange", y = moduleGenes["orange"], label = "*", color = "black", size = 8, vjust = 0)
ngenes
dev.off()

