# tpm function

tpm2 <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

tpm <-
  function(counts, len) {
    rpk <- counts / len
    return(t(t(rpk) * 1e6 / colSums(rpk)))
  }


# Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
# Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
# Divide the RPK values by the “per million” scaling factor. This gives you TPM.


# All genes
sub_list <- select(toptab, Description)
sub_list$hgnc_symbol <- sub_list$Description

# Top DGEs
HB <- c("HBA1", "HBA2", "HBB", "HBG1", "HBG2")
sub_list <- filter(toptab, (adj.P.Val < 0.05 & (logFC < -0.1 | logFC > 0.1 ))| Description %in% HB) %>% select(Description) # 299 - 162
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
top_pheno <- cbind(topmed$samples, top_tpm)
top_pheno <- select(top_pheno, -commonid)

# Corr matrix
control_vec <- top_tpm$hbb # pick HBB gene to correlate TPM values with
top_tpm <- as.matrix(top_tpm)

# calculate significant correlations
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





