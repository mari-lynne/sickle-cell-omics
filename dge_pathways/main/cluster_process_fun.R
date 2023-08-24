
# Function to process significant colors and generate gene_interest and overlap tables
process_color_module <- function(module_colors, gene_matrix, ensembl_mart, results_dir, top_table = NULL) {
  for (mod_col in module_colors) {
    # Get gene_interest table
    gene_interest <- colnames(gene_matrix)[module_colors == mod_col]
    gene_interest <- sub("\\.[0-9]+$", "", gene_interest)
    
    # Get gene names from Ensembl using the provided mart object
    genes <- getBM(
      attributes = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name'),
      filters = 'ensembl_gene_id',
      values = gene_interest,
      mart = ensembl_mart
    )
    
    # Write gene_interest table to CSV
    write.csv(genes, file = paste0(results_dir, "/gene_interest_", mod_col, ".csv"))
    
    # Process top_table if provided
    if (!is.null(top_table)) {
      genes_df <- as.data.frame(genes)
      toptab <- read.csv(top_table)
      toptab <- filter(toptab, (adj.P.Val < 0.05 & (logFC < -0.1 | logFC > 0.1 ))) %>%
        select(Description)
      
      # Get overlap table
      overlap <- genes_df[genes_df$hgnc_symbol %in% toptab$Description, ]
      write.csv(overlap, file = paste0(results_dir, "/overlap_", mod_col, ".csv"))
    }
  }
}
