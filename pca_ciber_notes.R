# TODO:

# Match sct pc's samples with rnaseq ids
# Run cibersort local script (in CSeQTL) for all sct RNAseq participants
# Add PCs and cell types into dge modelling

# PC's 

#/fh/scratch/delete90/kooperberg_c/sct_wgs/gds.minDP0/PCA
#scp mjohnso5@rhino1:/fh/scratch/delete90/kooperberg_c/sct_wgs/gds.minDP0/PCA/PC1_PC2.jpeg ~/Documents/CSeQTL/data

#  Cibersort

Check any latent batch effect in RNA-seq data.

Use log-transformed TReC, e.g, log(TReC of gene j of sample i/read-depth of sample i) as response variable,
fit a linear model with all the covariates included for eQTL mapping.
Then we need to look at the top a few PCs, check the variance they explain to decide how many PCs to use.
One way to check whether we need one covariate is to run the baseline model (the model without SNP genotype) for all genes,
and then check the histogram of the p-values across all genes for each covariate. If the histogram shows enrichment of small p-value,
then we should keep the covariate.
If the histogram shows the p-value distribution is uniform, then we may not need that covariate. 