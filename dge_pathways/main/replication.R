# Replication genes

name <- "/top_tab.csv"
data <- read.csv(file = paste0(wd, name)) # results_dir
top_data <- filter(data, adj.P.Val < 0.05)
trans2 <- filter(data, adj.P.Val < 0.05)

# SCA steady state genes 
sca_genes <- read.csv(file = "~/Documents/whi_sca/rna/sca_dge.csv")
sca_genes$Description <- str_trim(sca_genes$Description)
sig_sca <- filter(top_data, Description %in% sca_genes$Description)

# DGE SCA Crisis vs Healthy controls
sca_crisis <- read.csv(file = "~/Documents/whi_sca/rna/sca_crisis_dge.csv")
sca_crisis$Description <- str_trim(sca_crisis$Description)
sig_crisis <- filter(top_data, Description %in% sca_crisis$Description)

rm <- intersect(sig_crisis$Description, sig_sca$Description)
sca_all <- filter(sig_crisis, Description %!in% rm) %>% bind_rows(sig_sca)

# Enrichment:
test <- sig_sca
test$hgnc_symbol <- test$Description
genes <- test$hgnc_symbol

# Map gene symbol to enterez ids
genes_ez <- Organism.dplyr::select(src, 
                                   keys = genes,
                                   columns = c("symbol","entrez"),
                                   keytype = "symbol")
genes_ez <- distinct(genes_ez)

# add back to orignal df
genes_ez <- rename(genes_ez, hgnc_symbol = symbol)
test <- left_join(test, genes_ez, by = "hgnc_symbol")

ego <- enrichGO(gene = genes_ez$entrez,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", readable = TRUE)


View(ego@result)
dotplot(ego)



mesa <- fread("~/Documents/whi_sca/rna/MESA_qtl.csv")
mesa <- filter(mesa, pvalue <0.05) # 4335
mesa_whi_trans2 <- intersect(mesa$Gene_stable_ID, trans2$Name)
mesa_whi_trans2 <- trans2[trans2$Name %in% mesa$Gene_stable_ID, ]

# 363 relaxed filtering

mesa_whi_trans <- intersect(mesa$Gene_stable_ID, trans_sig$ensembl_gene_id)
mesa_whi_trans <- trans_sig[trans_sig$ensembl_gene_id %in% mesa_whi_trans, ]