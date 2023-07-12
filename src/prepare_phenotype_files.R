
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
  
  stop(paste0("Three command line arguments expected: regenie results file path. ", length(args), " found."))
  
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

delivery_table <- read.table(
  file = file.path(tables_folder, "delivery.gz"),
  header = T,
  sep = "\t"
) %>% 
  filter(
    !is.na(child_sentrix_id) | !is.na(mother_sentrix_id) | !is.na(father_sentrix_id)
  )

anthropometrics_table <- read.table(
  file = file.path(tables_folder, "child_anthropometrics.gz"),
  header = T,
  sep = "\t"
) %>% 
  filter(
    !is.na(child_sentrix_id) | !is.na(mother_sentrix_id) | !is.na(father_sentrix_id)
  ) %>% 
  select(
    -preg_id, -rank_siblings, -mother_id, -father_id, 
    -child_sentrix_id, -mother_sentrix_id, -father_sentrix_id, 
    -unrelated_children, -child_batch, -mother_batch, -father_batch, 
    -pregnancy_duration_term, -pregnancy_duration_preterm, -sex
  )

pregnancy_table <- read.table(
  file = file.path(tables_folder, "pregnancy.gz"),
  header = T,
  sep = "\t"
) %>% 
  filter(
    !is.na(child_sentrix_id) | !is.na(mother_sentrix_id) | !is.na(father_sentrix_id)
  ) %>% 
  select(
    -preg_id, -rank_siblings, -mother_id, -father_id, 
    -child_sentrix_id, -mother_sentrix_id, -father_sentrix_id, 
    -unrelated_children, -child_batch, -mother_batch, -father_batch, 
    -pregnancy_duration_term, -pregnancy_duration_preterm, -sex
  )

print(paste0(Sys.time(), " - Merging and exclusion criteria"))

pheno_table <- delivery_table %>% 
  inner_join(
    anthropometrics_table,
    by = "child_id"
  ) %>% 
  inner_join(
    pregnancy_table,
    by = "child_id"
  ) %>% 
  filter(
    pregnancy_duration_term == 1 & !is.na(umbilical_cord_length) & umbilical_cord_length > 0
  )

# Correct/Remove outliers

print(paste0(Sys.time(), " - Handling of outliers"))

pheno_table <- pheno_table %>% 
  mutate(
    umbilical_cord_length = ifelse(umbilical_cord_length >= 250, umbilical_cord_length / 10, umbilical_cord_length),
    umbilical_cord_length = ifelse(umbilical_cord_length <= 10, umbilical_cord_length * 10, umbilical_cord_length)
  )

mean <- mean(pheno_table$umbilical_cord_length)
sd <- sd(pheno_table$umbilical_cord_length)

pheno_table <- pheno_table %>% 
  mutate(
    umbilical_cord_length = ifelse(abs(umbilical_cord_length - mean) > 5 * sd, NA, umbilical_cord_length)
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
        sex,
        umbilical_cord_length
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
        sex,
        umbilical_cord_length
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
        sex,
        umbilical_cord_length
      ),
    by = "IID",
    multiple = "any"
  )

for (batch_name in unique(pheno_table_gwas_father$father_batch)) {
  
  batch_column <- make_clean_names(paste0("batch_", tolower(batch_name)))
  
  batch_column <- 
  
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



