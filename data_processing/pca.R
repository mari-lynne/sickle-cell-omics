# 23/03/23
# Aims
# Run PCA of expression data to account for population structure/variance in datra when modelling eQTLs
# For genome wide associations, it may be better to include genetic PCs 

library("FactoMineR")
library("factoextra")
library("PCAForQTL")

setwd("~/Documents/whi_sca/rna/results")
# load(file = "prelim_dge.RData")

exprs <- t(dge$counts)

prcompResult <-prcomp(exprs,center=TRUE,scale.=TRUE)
PCs <-prcompResult$x
dim(PCs)

resultRunElbow <-PCAForQTL::runElbow(prcompResult=prcompResult)
print(resultRunElbow)


K_elbow <-resultRunElbow #12.
K_BE <-resultRunBE$numOfPCsChosen #29.
K_GTEx <-60 #GTEx uses 60 PEER factors, and they are almost identical to the top 60 PCs.
PCAForQTL::makeScreePlot(prcompResult,labels=c("Elbow"),values=c(K_elbow),
                         titleText="RNA-Seq PCA")
ggsave(filename = "scree.png", device = "png")
pc_top <-PCs[,1:K_elbow] 

write.table(pc_top, file = "pc_top.txt", quote = F, col.names = T, row.names = T)

# update covar matrix (do at the start of sca_deg.R)

# check against known covars ---------------------------


knownCovariates<- select(dge$samples, "age", "bmi_t0")
  

dataCovariates[,c(1:5,66:68)] #368*8. 368 samples, 8 known covariates.
identical(rownames(knownCovariates),rownames(expr))


covars <- cbind(pheno, pc_top)
