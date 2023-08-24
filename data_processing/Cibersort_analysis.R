Deconvolution -------------------------------------------------------------

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

# Merge cell types -------------------------------------------------------------
pp_hat_ciber <- read.csv(file = "~/Documents/whi_sca/rna/results/ciber/R/pp_hat_ciber_test.csv")

# Merged column names:
T_cells <- colnames(pp_hat_ciber)[str_starts(colnames(pp_hat_ciber), "T")==TRUE]
B_cells <- colnames(pp_hat_ciber)[str_starts(colnames(pp_hat_ciber), "B|Plasma")==TRUE]
NK <- colnames(pp_hat_ciber)[str_starts(colnames(pp_hat_ciber), "NK")==TRUE]
mono <- colnames(pp_hat_ciber)[str_starts(colnames(pp_hat_ciber), "Mono|Macro|Den|Mast")==TRUE]
neut <- colnames(pp_hat_ciber)[str_starts(colnames(pp_hat_ciber), "Neut|Eo")==TRUE]

# Create a new matrix with combined values and columns
combined_matrix <- matrix(0, nrow = nrow(pp_hat_ciber), ncol = 5)  # Adjust the number of columns as needed
colnames(combined_matrix) <- c("T_cells", "B_cells", "NK", "mono", "neut")
row.names(combined_matrix) <- row.names(pp_hat_ciber)
row.names(combined_matrix) <- (pp_hat_ciber$X)

combined_matrix[, "T_cells"] <- rowSums(pp_hat_ciber[, T_cells])
combined_matrix[, "B_cells"] <- rowSums(pp_hat_ciber[, B_cells])
combined_matrix[, "NK"] <- rowSums(pp_hat_ciber[, NK])
combined_matrix[, "mono"] <- rowSums(pp_hat_ciber[, mono])
combined_matrix[, "neut"] <- rowSums(pp_hat_ciber[, neut])

write.csv(combined_matrix, file = "~/Documents/whi_sca/rna/results/ciber/R/pp_hat_ciber_merged.csv")


# Plotting ---------------------------------------------------------------------

# ciber <- as.data.frame(pp_hat_ciber)
ciber <- as.data.frame(combined_matrix)
ciber$sample_id <- row.names(ciber)

ciber <- ciber %>% group_by(sample_id)
ciber <- melt(ciber, id.vars = c("sample_id"), measure.vars = c(1:(ncol(ciber)-1)))
ciber <- ciber %>% rename(fraction = value,
                          cell_type = variable)
# Plot
ciber %>% ggplot(aes(x=sample_id, y=fraction, fill=cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete=TRUE, option = "H", direction = -1, alpha = 0.95) +
  theme_light() +
  labs(y = "Fraction \n", x= "Sample ID", title = "Cibersort Output")


