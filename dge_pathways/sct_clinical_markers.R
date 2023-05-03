# Plan

# Use top dge genes and see if they stratify with egfr/cell counts or can predict it

# Kidney data

whi <- read.csv(file = paste0(meta_dir, "/WHI_HTNKidney_2023-02-27.csv"))
whi <- rename(whi, subject_id = subjectID)

tal1 <- inner_join(tal1, whi, by = "subject_id")

# TAL1

tal1_sct <- filter(tal1, sct == "1")
table(tal1_sct$ckd) # 12 ckd
table(tal1_sct$eGFRMDRD)

tal1_sct %>% ggplot(aes(x = eGFRCKDEpi, y = TAL1_counts)) +
  geom_point() +
  theme_light() +
  scale_y_continuous(trans='log10')



tal1 %>% filter(sct == "1") %>%
  ggplot(aes(x = sct, y = TAL1_counts)) +
  geom_boxplot()+
  theme_light() +
  scale_y_continuous(trans='log10') +
  labs(y = "TAL1 expression (counts)\n", x = "SCT status") +
  stat_compare_means(comparisons = comp)

