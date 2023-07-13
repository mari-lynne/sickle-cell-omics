# 5) Deconvolution -------------------------------------------------------------

# AIMS:
# Using trecase output from step 4; deconvolute cell fractions using cibersort and LM22 signature matrix

## Dirs and file names ----------------------------------------------------------

library(data.table)
library(dplyr)
library(stringr)
library(smarter)

work_dir <- c("~/Documents/whi_sca/rna") # Where main data is saved
plot_dir <- c("~/Documents/whi_sca/rna/plots") # Where to save output plots
id_dir <- ("~/Documents/whi_sca/rna/ids")
meta_dir <- c("~/Documents/whi_sca/rna/meta")
results_dir <- c("~/Documents/whi_sca/rna/results")

setwd(work_dir)

# Input files
sig_tpm <- fread("~/Documents/CSeQTL/data/ciber_ase/LM22.txt")
bulk <- fread(paste0(results_dir, "/sct_tpm_rnaseq.txt"))

# Output files
sig_fn = file.path(results_dir, "signature.txt") # TPM signature expression matrix
mix_fn = file.path(results_dir, "mixture.txt") # TPM bulk expression matrix

# write tables
write.table(sig_tpm, file = sig_fn, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(bulk, file = mix_fn, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

## Run cibersort
source("~/Documents/CSeQTL/scripts/CSeQTL/5_cibersort_source.R") # Obtained from CIBERSORT website

work_dir = "~/Documents/CSeQTL/data/ciber_ase"
sig_fn = file.path(work_dir, "signature.txt") # TPM signature expression matrix
mix_fn = file.path(work_dir, "mixture.txt") # TPM bulk expression matrix

# Output 
results <- CIBERSORT(sig_matrix = sig_fn, mixture_file = mix_fn, filename = "DECON",
                     perm = 0, QN = FALSE, absolute = FALSE, abs_method = 'sig.score')

# filename = "DECON"
results = as.data.frame(results)
# Extract proportion of transcripts per cell type per sample
pp_bar_ciber <- results[1:22] # 22 immune subsets

# Cell size table (estimated in EPIC paper) 
# TODO update with our own FSC/SSC data
T_cells <- colnames(pp_bar_ciber)[str_starts(colnames(pp_bar_ciber), "T")==TRUE]
B_cells <- colnames(pp_bar_ciber)[str_starts(colnames(pp_bar_ciber), "B|Plasma")==TRUE]
NK <- colnames(pp_bar_ciber)[str_starts(colnames(pp_bar_ciber), "NK")==TRUE]
mono <- colnames(pp_bar_ciber)[str_starts(colnames(pp_bar_ciber), "Mono|Macro|Den|Mast")==TRUE]
neut <- colnames(pp_bar_ciber)[str_starts(colnames(pp_bar_ciber), "Neut|Eo")==TRUE]

sizes <- c(rep(0.4, length(B_cells)),
               rep(0.4, length(T_cells)),
               rep(0.42, length(NK)),
               rep(1.4, length(mono)),
               rep(0.15, length(neut)))

pp_bar_ciber <- rbind(pp_bar_ciber, sizes)

# Adjust expression for cell size S
pp_hat_ciber <- t(apply(pp_bar_ciber, 1, function(xx) {
  yy <- xx / sizes
  yy / sum(yy)
}))

pp_hat_ciber <- pp_hat_ciber[-nrow(pp_hat_ciber),]

write.csv(pp_hat_ciber, file = paste0(results_dir, "pp_hat_ciber_test.csv"))

# Plotting ---------------------------------------------------------------------

ciber <- as.data.frame(pp_hat_ciber)
# ciber <- ciber %>% mutate(`941526` = jitter(`941437`))

test <- as.data.frame(pp_hat_ciber)
test$sample_id <- row.names(test)

test <- test %>% group_by(sample_id)
test <- melt(test, id.vars = c("sample_id"), measure.vars = c(1:(ncol(test)-1)))
test <- test %>% rename(fraction = value,
                        cell_type = variable)

library(viridis)
library(ggplot2)
test %>% ggplot(aes(x=sample_id, y=fraction, fill=cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete=TRUE, option = "H", direction = -1, alpha = 0.95) +
  theme_light() +
  labs(y = "Fraction \n", x= "Sample ID", title = "Cibersort Output Test") %>% facet_wrap(~)

ciber %>% ggplot(aes())
