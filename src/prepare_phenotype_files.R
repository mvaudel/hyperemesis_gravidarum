
##
#
# This script prepares phenotype files for Regenie
#
##

# Seed for random choice of sample

set.seed(05072023)

# Libraries

library(conflicted)
library(janitor)
library(glue)
library(dplyr)
library(ggplot2)
library(grid)

conflicts_prefer(dplyr::filter)


# General parameters
theme_set(theme_bw(base_size = 16))


# Command line input

args <- commandArgs(TRUE)

if (length(args) != 3) {
  
  stop(paste0("Three command line arguments expected. ", length(args), " found."))
  
}

tables_folder <- args[1]

if (!dir.exists(tables_folder)) {
  
  stop(paste0("Phenotype tables folder ", tables_folder, " not found."))
  
}

gwas_pheno_folder <- args[2]

if (!dir.exists(gwas_pheno_folder)) {
  
  stop(paste0("GWAS phenotypes folder ", gwas_pheno_folder, " not found."))
  
}

fam_file <- args[3]

if (!file.exists(fam_file)) {
  
  stop(paste0("Fam file ", fam_file, " not found."))
  
}


# Pheno files

print(paste0(Sys.time(), " - Loading phenotypes from '", tables_folder, "'"))

pregnancy_table <- read.table(
  file = "/mnt/archive/moba/pheno/v12/pheno_anthropometrics_23-07-12/pregnancy.gz",
  header = T,
  sep = "\t"
)

# Exclusion criteria

print(paste0(Sys.time(), " - Exclusion criteria"))

pregnancy_table <- pregnancy_table %>%
  filter(
    is.na(plural_birth) & pregnancy_duration > 90
  )

# Format variables of interest

print(paste0(Sys.time(), " - Hadnling of outliers and combination of variables"))

pheno_table <- pregnancy_table %>% 
  mutate(
    hospitalized_prolonged_nausea_vomiting = ifelse(is.na(hospitalized_prolonged_nausea_vomiting), 0, 1),
    nausea_q2 = ifelse(!is.na(nausea_q2) & nausea_q2 == "Yes", 1, 0),
    vomiting_q2 = ifelse(!is.na(vomiting_q2) & vomiting_q2 == "Yes", 1, 0),
    nausea_any = ifelse(nausea_q2 == 1 | !is.na(nausea_before_4w) | !is.na(nausea_5w_8w) | !is.na(nausea_9w_12w) | !is.na(nausea_13w_15w) | !is.na(nausea_13w_16w) | !is.na(nausea_17w_20w) | !is.na(nausea_21w_24w) | !is.na(nausea_25w_28w) | !is.na(nausea_after_29w), 1, 0),
    vomiting_any = ifelse(vomiting_q2 == 1 | !is.na(vomiting_before_4w) | !is.na(vomiting_5w_8w) | !is.na(vomiting_9w_12w) | !is.na(vomiting_13w_15w), 1, 0),
    long_term_nausea_vomiting_any = ifelse(!is.na(long_term_nausea_vomiting_13w_16w) | !is.na(long_term_nausea_vomiting_17w_20w) | !is.na(long_term_nausea_vomiting_21w_24w) | !is.na(long_term_nausea_vomiting_25w_28w) | !is.na(long_term_nausea_vomiting_after_29w), 1, 0),
    nausea_vomiting = ifelse(nausea_any == 0 & vomiting_any == 0 & long_term_nausea_vomiting_any == 0, 0, 1),
    nausea_duration = ifelse(is.na(nausea_week_most_bothered_from_q2) | is.na(nausea_week_most_bothered_to_q2), NA, nausea_week_most_bothered_to_q2 - nausea_week_most_bothered_from_q2),
    vomiting_duration = ifelse(is.na(vomiting_week_from_q2) | is.na(vomiting_week_to_q2), NA, vomiting_week_to_q2 - vomiting_week_from_q2),
    nausea_week_most_bothered_from_q2 = ifelse(!is.na(nausea_duration) & nausea_duration < 0, NA, nausea_week_most_bothered_from_q2),
    nausea_week_most_bothered_to_q2 = ifelse(!is.na(nausea_duration) & nausea_duration < 0, NA, nausea_week_most_bothered_to_q2),
    nausea_duration = ifelse(!is.na(nausea_duration) & nausea_duration < 0, NA, nausea_duration),
    vomiting_week_from_q2 = ifelse(!is.na(vomiting_duration) & vomiting_duration < 0, NA, vomiting_week_from_q2),
    vomiting_week_to_q2 = ifelse(!is.na(vomiting_duration) & vomiting_duration < 0, NA, vomiting_week_to_q2),
    vomiting_duration = ifelse(!is.na(vomiting_duration) & vomiting_duration < 0, NA, vomiting_duration)
  ) %>% 
  mutate(
    hg_vs_all = ifelse(hospitalized_prolonged_nausea_vomiting == 1, 1, 0),
    hg_vs_no_nausea_vomiting = ifelse(nausea_vomiting == 0, 0, NA),
    hg_vs_no_nausea_vomiting = ifelse(hospitalized_prolonged_nausea_vomiting == 1, 1, hg_vs_no_nausea_vomiting),
    vomiting_before_4w = ifelse(!is.na(vomiting_before_4w), 1, 0),
    vomiting_5w_8w = ifelse(!is.na(vomiting_before_4w), 1, 0),
    vomiting_9w_12w = ifelse(!is.na(vomiting_9w_12w), 1, 0),
    vomiting_13w_15w = ifelse(!is.na(vomiting_13w_15w), 1, 0),
    long_term_nausea_vomiting_13w_16w = ifelse(!is.na(long_term_nausea_vomiting_13w_16w), 1, 0),
    long_term_nausea_vomiting_17w_20w = ifelse(!is.na(long_term_nausea_vomiting_17w_20w), 1, 0),
    long_term_nausea_vomiting_21w_24w = ifelse(!is.na(long_term_nausea_vomiting_21w_24w), 1, 0),
    long_term_nausea_vomiting_25w_28w = ifelse(!is.na(long_term_nausea_vomiting_25w_28w), 1, 0),
    long_term_nausea_vomiting_after_29w = ifelse(!is.na(long_term_nausea_vomiting_after_29w), 1, 0)
  )


# Fam file

print(paste0(Sys.time(), " - Loading fam file"))

fam_table <- read.table(
  file = fam_file,
  header = F,
  sep = " "
)
fam_table <- fam_table[, 1:2]
names(fam_table) <- c("FID", "IID")


# Set up data frames for GWAS

print(paste0(Sys.time(), " - Setting up for gwas"))

pheno_table_gwas_child <- fam_table %>% 
  inner_join(
    pheno_table %>%
      filter(
        !is.na(child_sentrix_id)
      ) %>% 
      select(
        IID = child_sentrix_id,
        child_batch,
        hg_vs_all,
        hg_vs_no_nausea_vomiting,
        vomiting_before_4w,
        vomiting_5w_8w,
        vomiting_9w_12w,
        vomiting_13w_15w,
        long_term_nausea_vomiting_13w_16w,
        long_term_nausea_vomiting_17w_20w,
        long_term_nausea_vomiting_21w_24w,
        long_term_nausea_vomiting_25w_28w,
        long_term_nausea_vomiting_after_29w,
        vomiting_week_from_q2,
        vomiting_week_to_q2,
        vomiting_duration
      ),
    by = "IID"
  )

for (batch_name in unique(pheno_table_gwas_child$child_batch)) {
  
  batch_column <- make_clean_names(paste0("batch_", tolower(batch_name)))
  
  pheno_table_gwas_child[[batch_column]] <- ifelse(pheno_table_gwas_child$child_batch == batch_column, 1, 0)
  
}

pheno_table_gwas_mother <- fam_table %>% 
  inner_join(
    pheno_table %>% 
      filter(
        !is.na(mother_sentrix_id)
      ) %>% 
      select(
        IID = mother_sentrix_id,
        mother_batch,
        hg_vs_all,
        hg_vs_no_nausea_vomiting,
        vomiting_before_4w,
        vomiting_5w_8w,
        vomiting_9w_12w,
        vomiting_13w_15w,
        long_term_nausea_vomiting_13w_16w,
        long_term_nausea_vomiting_17w_20w,
        long_term_nausea_vomiting_21w_24w,
        long_term_nausea_vomiting_25w_28w,
        long_term_nausea_vomiting_after_29w,
        vomiting_week_from_q2,
        vomiting_week_to_q2,
        vomiting_duration
      ),
    by = "IID",
    multiple = "any"
  )

for (batch_name in unique(pheno_table_gwas_mother$mother_batch)) {
  
  batch_column <- make_clean_names(paste0("batch_", tolower(batch_name)))
  
  pheno_table_gwas_mother[[batch_column]] <- ifelse(pheno_table_gwas_mother$mother_batch == batch_column, 1, 0)
  
}

pheno_table_gwas_father <- fam_table %>% 
  inner_join(
    pheno_table %>% 
      filter(
        !is.na(father_sentrix_id)
      ) %>% 
      select(
        IID = father_sentrix_id,
        father_batch,
        hg_vs_all,
        hg_vs_no_nausea_vomiting,
        vomiting_before_4w,
        vomiting_5w_8w,
        vomiting_9w_12w,
        vomiting_13w_15w,
        long_term_nausea_vomiting_13w_16w,
        long_term_nausea_vomiting_17w_20w,
        long_term_nausea_vomiting_21w_24w,
        long_term_nausea_vomiting_25w_28w,
        long_term_nausea_vomiting_after_29w,
        vomiting_week_from_q2,
        vomiting_week_to_q2,
        vomiting_duration
      ),
    by = "IID",
    multiple = "any"
  )

for (batch_name in unique(pheno_table_gwas_father$father_batch)) {
  
  batch_column <- make_clean_names(paste0("batch_", tolower(batch_name)))
  
  pheno_table_gwas_father[[batch_column]] <- ifelse(pheno_table_gwas_father$father_batch == batch_column, 1, 0)
  
}


# Write tables

print(paste0(Sys.time(), " - Export of tables"))

write.table(
  x = pheno_table_gwas_child,
  file = file.path(gwas_pheno_folder, "pheno_child"),
  row.names = F,
  col.names = T,
  quote = F
)

write.table(
  x = pheno_table_gwas_mother,
  file = file.path(gwas_pheno_folder, "pheno_mother"),
  row.names = F,
  col.names = T,
  quote = F
)

write.table(
  x = pheno_table_gwas_father,
  file = file.path(gwas_pheno_folder, "pheno_father"),
  row.names = F,
  col.names = T,
  quote = F
)



