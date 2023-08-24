setwd("/home/mari/Documents/whi_sca")

jeff <- read.csv(file = "whi_sct_wgs_id_map.csv")
index <- read.csv(file = "wgs_cram_index.csv")

check <- intersect(index$file_prefix, jeff$SAMPLE_ID) # 349 samples
sct_overlap <- jeff[jeff$SAMPLE_ID %in% check,]

meta <- read.csv(file = "~/Documents/whi_sca/rna/meta/sct_all_covars_Jul23.csv")
check <- intersect(sct_overlap$SUBJECT_ID, meta$subject_id) # 278

sct_meta <- meta[meta$subject_id %in% check,]
