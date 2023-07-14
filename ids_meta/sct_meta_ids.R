### Aims:

# Gather all IDs and phenotype information for SCT rna-seq data
# Add clinical data
# Produce 'clean' covariate files for downstream analysis

### R set up and Functions -----------------------------------------------------

source("~/scripts/r-packages.R")
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

# Match RNAseq IDs to nwgcid_topmed IDs

# Original id and link files copied from /fh/scratch/delete90/kooperberg_c/sct_rnaseq
# and /fh/fast/kooperberg_c/OMICS/id_maps/WHI_GWAS-WGS_OMICS_overlap_2022-12-16.csv
# moved to local folder ~/Documents/whi_sca/rna/ids

id_link <- clean_names(read.csv(file = "~/Documents/whi_sca/rna/ids/lookup_whi_topmed_to6_rnaseq_1_final.csv"))

# RNA-seq sample ids
rna_ids <- as.data.frame(colnames(data))
# Remove X string
rna_ids$nwgc_sample_id <- as.integer(str_remove(rna_ids$`colnames(data)`, "X"))
rna_link <- inner_join(rna_ids, id_link, by = "nwgc_sample_id")

# Topmed link file
id_link2 <- clean_names(read.csv(file = paste(id_dir, "as666_sct_rnaseq_ids.csv", sep = "/")))
rna_link <- inner_join(rna_link, id_link2, by = "torid") # excludes KS Sentinel samples
str(rna_link)

#### Tidy variables -------------------------------------------------------------

# rename plate vars

rna_link <- rename(rna_link, pick_plate = pick_plate_1_indicates_which_samples_were_prepped_together)

rna_link <- rna_link %>%
  mutate(plate = coalesce(rna_link$pick_plate, rna_link$pick_plate_2_repeat)) %>%
  select(-c("pick_plate", "pick_plate_2_repeat"))

rna_link <- rename(rna_link, flagged_qc = flagged_seq_qc_metric)
rna_link <- rename(rna_link, rnaseq_ids = `colnames(data)`)

# Remove flagged QC samples
rna_link <- filter(rna_link, flagged_qc != "Yes")

# 3) Updated LLS phenotype files -----------------------------------------------

# See Jeff correspondance
# ID to phenotype link file

id_link3 <- clean_names(read.csv(file = paste(id_dir, "WHI_GWAS-WGS_OMICS_overlap_2022-12-16.csv", sep = "/")))
rna_link <- inner_join(rna_link, id_link3, by = "subject_id")


#### RNAseq PCs ----------------------------------------------------------------

# See pca.R
pc_top <- fread(file = paste0(results_dir, "/pc_top.txt"))
pc_top <- rename(pc_top, rnaseq_ids = V1)
covar <- left_join(rna_link, pc_top, by = "rnaseq_ids")


# 4) clinical phenotype files ---------------------------------------------------

# Add whi data
whi <-
  read.csv(file = paste0(meta_dir, "/WHI_HTNKidney_2023-02-27.csv")) %>% 
  rename(subject_id = subjectID)

# Join by subject id
whi_sct <- filter(whi, subject_id %in% covar$subject_id)
covar <- left_join(covar, whi_sct, by= "subject_id") %>% select(where(not_all_na))
covar <- rename(covar, eGFR = eGFRCKDEpi)

# BMI phenotypes etc. 
whi_sct <- clean_names(fread(file = paste0(meta_dir, "/WHI_all_subjects_PAGE_updated_phenotypes.txt"))) %>% rename(subject_id = subjectid)

covar <- left_join(covar, whi_sct, by = "subject_id") # matched rows 837

# Clinical phenotypes + cell counts 
whi_clin <- clean_names(read.csv(file = paste0(meta_dir, "/WHI_inflamation_2023-02-27.csv")))
whi_clin <- filter(whi_clin, subject_id %in% covar$subject_id)

covar <- left_join(covar, whi_clin, by = "subject_id")

# SCT status -----------------------------------------------------------------
sct <- clean_names(read.csv(file = "~/Documents/whi_sca/rna/meta/sct_trait_info.csv"))

# Join and clean all duplicate columns
covar <- clean_cols(left_join(covar, sct, by = "commonid"))
# 192 no info, 837 sct info (no info are probably the hispanics, exclude for now)
covar$sct <- ifelse(is.na(covar$sct) == T, 0, 1)
table(covar$sct)

save.image(paste0(results_dir, "/metadata_july23.RData"))
write.csv(covar, file = paste0(meta_dir, "/sct_all_covars_Jul23.csv"), row.names = F)

# HBB ratios -------------------------------------------------------------------

# See 1_rnaseq_preprocess.R


# 5) Make final phenotype files and select vars --------------------------------

sct_pheno2 <- sct_pheno2 %>% select(rnaseq_ids, commonid, subject_id, bmi_t0, smoking, smokestatus, age.y, ethnic, sct, whills, flagged_seq_qc, pick_plate1, pick_plate_2_repeat)
sct_pheno2 <- sct_pheno2 %>% rename(age = age.y)

write.table(sct_pheno2, file = paste0(meta_dir, "/sct_rnaseq_pheno.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
write.table(sct_pheno2, file = paste0(meta_dir, "/rnaseq_pheno_all.txt"), quote = F, sep = "\t", row.names = F, col.names = T)

