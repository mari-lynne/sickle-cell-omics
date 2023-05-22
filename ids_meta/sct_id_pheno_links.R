### Aims:

# Gather all IDs and phenotype information for SCT rna-seq data

### R set up and Functions -----------------------------------------------------

# Bioconductor packages
library(BiocManager)
library(edgeR)
library(limma)

# Data cleaning
library(readODS)
library(forcats)
library(stringi)
library(stringr)
library(janitor)
library(data.table)
library(dplyr)
library(tidylog)

source("~/scripts/functions.R")

# Options ---------------------------------------------------------------------

wd <- c("~/Documents/whi_sca/rna") # Where main data is saved
plot_dir <- c("~/Documents/whi_sca/rna/plots") # Where to save output plots
id_dir <- ("~/Documents/whi_sca/rna/ids")
meta_dir <- c("~/Documents/whi_sca/rna/meta")
results_dir <- c("~/Documents/whi_sca/rna/results")

data_name <- c("whi_topmed_to6_rnaseq_gene_reads.gct.gz") 
pheno_name <- c("old/rna_seq_prelim_phenotype_2020-03-05.ods")

# 1) Data Set up ---------------------------------------------------------------

# RNA seq
data <- read_omic(name = data_name, wd = wd)

# 2) ID matching ---------------------------------------------------------------

# Orginal id and link files copied from /fh/scratch/delete90/kooperberg_c/sct_rnaseq
# and /fh/fast/kooperberg_c/OMICS/id_maps/WHI_GWAS-WGS_OMICS_overlap_2022-12-16.csv
# moved to local folder ~/Documents/whi_sca/rna/ids

id_link <- clean_names(read.csv(file = "~/Documents/whi_sca/rna/ids/lookup_whi_topmed_to6_rnaseq_1_final.csv"))
# id_link$subjectID <- as.character(id_link$subjectID)

# RNA-seq sample ids
rna_ids <- as.data.frame(colnames(data))
# Remove X string
rna_ids$nwgc_sample_id <- as.integer(str_remove(rna_ids$`colnames(data)`, "X"))
rna_link <- inner_join(rna_ids, id_link, by = "nwgc_sample_id")

# Topmed link file
id_link2 <- clean_names(read.csv(file = paste(id_dir, "as666_sct_rnaseq_ids.csv", sep = "/")))
rna_link <- inner_join(rna_link, id_link2, by = "torid") # excludes KS Sentinel samples

# ID to phenotype link file
id_link3 <- clean_names(read.csv(file = paste(id_dir, "WHI_GWAS-WGS_OMICS_overlap_2022-12-16.csv", sep = "/")))
rna_link <- inner_join(rna_link, id_link3, by = "subject_id")

# 3) Phenotype info ------------------------------------------------------------
# Downloaded from Jeff/Alex email chain
sct <- read.csv(file = "~/Documents/whi_sca/rna/meta/sct_trait_info.csv")

check <- inner_join(rna_link, sct, by = "commonid")
# 192 no info, 837 sct info (no info are probably the hispanics, exclude for now)
check$sct <- ifelse(is.na(check$sct) == T, 0, 1)
table(check$sct)
sct_pheno <- rna_link #check

# ID list
samples <- rna_link$commonid
write.table(samples, file = paste(id_dir, "commonid_sct_rnaseq.txt", sep = "/"), quote = F, sep = "\t", row.names = F, col.names = F)

# Tidy pheno file
sct_pheno <- sct_pheno %>% dplyr::rename(rnaseq_ids = 'colnames(data)')
sct_pheno <- sct_pheno %>% dplyr::select(rnaseq_ids, commonid, subject_id, age, ethnic, flagged_seq_qc_metric, pick_plate_1_indicates_which_samples_were_prepped_together, pick_plate_2_repeat)
colnames(sct_pheno) <- str_remove_all(colnames(sct_pheno), "\\.y")
sct_pheno <- sct_pheno %>% rename(flagged_seq_qc = flagged_seq_qc_metric,
                                  pick_plate1 = pick_plate_1_indicates_which_samples_were_prepped_together)

write.table(sct_pheno, file = paste0(meta_dir, "/sct_rnaseq_pheno.txt"), quote = F, sep = "\t", row.names = F, col.names = T)

# 3b) Extra pheno info - BMI ---------------------------------------------------

meta <- clean_names(fread(file = paste0(meta_dir, "/WHI_all_subjects_PAGE_updated_phenotypes.txt")))
meta <- meta %>% rename(subject_id = subjectid)

sct_pheno2 <- left_join(sct_pheno, meta, by = "subject_id") # matched rows 837
sct_pheno2 <- sct_pheno2 %>% select(rnaseq_ids, commonid, subject_id, bmi_t0, smoking, smokestatus, age.y, ethnic, sct, whills, flagged_seq_qc, pick_plate1, pick_plate_2_repeat)
sct_pheno2 <- sct_pheno2 %>% rename(age = age.y)

write.table(sct_pheno2, file = paste0(meta_dir, "/sct_rnaseq_pheno.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
write.table(sct_pheno2, file = paste0(meta_dir, "/rnaseq_pheno_all.txt"), quote = F, sep = "\t", row.names = F, col.names = T)


# Manual genotype info attempts ------------------------------------------------
# Plink ids
# see bash_notes; ran plink code in terminal

# plink2 --vcf WHI_share_Affy6.0-2015-03-05.chr11.vcf.gz --make-pgen--out whi_share_chr11
psam <- fread("~/Documents/whi_sca/genomics/whi_share_chr11.psam")
keep <- intersect(psam$`#IID`, rna_link$share_sample_id)

write.table(keep, file = paste(id_dir, "keep_list.txt", sep = "/"), quote = F, sep = "\t", row.names = F, col.names = F)

# plink2 --pfile whi_share_chr11 --keep keep_list.txt --make-pgen --out whi_share_rna

# filter for sca variant # CHR 11 5227002

pvar <- fread("~/Documents/whi_sca/genomics/whi_share_rna.pvar")
sct <- filter(pvar, POS == "5227001") # rs311 site

sct <- filter(pvar, POS %in% 5226001:5228001)
# No variant measured at that site, possibly filtered in QC