# Set up ------

source("~/scripts/functions.R")
source("~/scripts/r-packages.R")

# Options and Directories ------------------------------------------------------

wd <- c("~/Documents/whi_sca/rna") # Where main data is saved
plot_dir <- c("~/Documents/whi_sca/rna/plots") # Where to save output plots
id_dir <- ("~/Documents/whi_sca/rna/ids")
meta_dir <- c("~/Documents/whi_sca/rna/meta")
results_dir <- c("~/Documents/whi_sca/rna/results")

# read in data 

setwd(wd)
load(file = paste0(results_dir, "/sct_dge.RData"))

# TPM --------------------------------------------------------------------------

# 1) Calculate TPM matrix ------------------------------------------------------

# Top DGEs
HB <- c("HBA1", "HBA2", "HBB", "HBG1", "HBG2")
# sub_list <- filter(toptab, (adj.P.Val < 0.05 & (logFC < -0.1 | logFC > 0.1 ))| Description %in% HB) %>% select(Description) # 299 - 162
sub_list <- select(toptab, Description)
sub_list$hgnc_symbol <- sub_list$Description

# Get gene/transcript length from biomart
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(
  attributes=c('ensembl_gene_id', 'hgnc_symbol','start_position','end_position'),
  mart = ensembl)
genes$size = genes$end_position - genes$start_position

# filter dups
sub_list <- sub_list[(stri_duplicated(sub_list$hgnc_symbol)==FALSE), ]
genes <- genes[(stri_duplicated(genes$hgnc_symbol)==FALSE), ]
# Join gene length info
sub_list <- filter(genes, hgnc_symbol %in% sub_list$hgnc_symbol)
dim(sub_list)

# Filter count data
topmed <- dge[dge$genes$Description %in% sub_list$hgnc_symbol, ]
topmed$genes$hgnc_symbol <- topmed$genes$Description
dim(topmed)
topmed$genes <- left_join(sub_list, topmed$genes, by = "hgnc_symbol")
dim(topmed$counts)

options(scipen=999)
top_tpm <- as.data.frame(t(tpm(topmed$counts, topmed$genes$size)))
colnames(top_tpm) <- topmed$genes$hgnc_symbol
top_tpm <- clean_names(top_tpm)

# Alex tables
# top_pheno <- cbind(topmed$samples, top_tpm)
# top_pheno <- select(top_pheno, -commonid)

# Corr matrix
control_vec <- top_tpm$hbb # pick HBB gene to correlate TPM values with
top_tpm <- as.matrix(top_tpm)

# Calculate all correlations --------------------------------------------------

correlations <- apply(top_tpm, 2, function(x) cor(x, control_vec,  method = "pearson"))
# Output correlations
print(correlations)

# Create an empty dataframe
results_df <- data.frame(Gene_ID = character(), Correlation = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)

# Compute correlation and p-values
for (i in 1:ncol(top_tpm)) {
  
  # subset data
  gene_expression <- top_tpm[,i ]
  # calculate correlations
  cor_result <- cor.test(gene_expression, control_vec, method = "pearson", adjust="none")
  correlation <- cor_result$estimate
  p_value <- cor_result$p.value
  
  # Add results to the dataframe
  results_df <- rbind(results_df, data.frame(Gene_ID = colnames(top_tpm)[i], Correlation = correlation, P_Value = p_value))
}

# Print the dataframe
print(results_df)
str(results_df)
results_sig <- filter(results_df, P_Value <= 0.05)
results_sig <- filter(results_df, P_Value <=0.05 & (Correlation > 0.5 | Correlation < -0.5))

all_cors <- results_df


correlations <- apply(top_tpm, 2, function(x) cor(x, control_vec,  method = "pearson"))
# Output correlations
print(correlations) 
correlations <- apply(top_tpm, 2, function(x) cor.test(x, control_vec,  method = "pearson"))


### Calculate significant correlations ------------------------------------

# Create an empty dataframe
results_df <- data.frame(Gene_ID = character(), Correlation = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)

# Compute correlation and p-values
for (i in 1:ncol(top_tpm)) {
  
  # subset data
  gene_expression <- top_tpm[,i ]
  # calculate correlations
  cor_result <- cor.test(gene_expression, control_vec, method = "pearson", adjust="none")
  correlation <- cor_result$estimate
  p_value <- cor_result$p.value
  
  # Add results to the dataframe
  results_df <- rbind(results_df, data.frame(Gene_ID = colnames(top_tpm)[i], Correlation = correlation, P_Value = p_value))
}

# Print the dataframe
print(results_df)
str(results_df)
results_sig <- filter(results_df, P_Value <= 0.05)
results_sig <- filter(results_df, P_Value <=0.05 & (Correlation > 0.5 | Correlation < -0.5))
all_cors <- results_df

### Plot sig correlations -------------------------------------------------

correlations <- apply(top_tpm, 2, function(x) {
  cor_result <- cor(x, control_vec, method = "pearson")
  
  # Filter correlations based on threshold and statistical significance
  if (abs(cor_result) > 0.5 && cor_result >= 0 && cor_result < 1) {
    p_value <- cor.test(x, control_vec, method = "pearson", exact = FALSE)$p.value
    if (p_value < 0.05) {
      # Plot correlation
      var_name <- names(cor_result)[x]
      df <- data.frame(x = log(control_vec), y = log(x))
      p <- ggplot(df, aes(x = x, y = y)) + 
        geom_smooth(method = "lm", fill = "#A92D2F", alpha = 0.7) +
        geom_point() +
        geom_text(aes(label = paste("r =", round(cor_result, 2)," p =", signif(p_value, 3))), x = Inf, y = -Inf, 
                  hjust = 2, vjust = 0, size = 4, color = "#d15224") + theme_minimal() +
        labs(x=control_vec_lab, y = var_name)
      print(p)
      
      return(cor_result)
    }
  }
  
  return(NA)  # If correlation doesn't meet the criteria, return NA
})

# Remove NA values
cor_result <- correlations[!is.na(correlations)]

# Output correlations
print(correlations)

