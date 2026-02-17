
##
#
# Merges results from child, mother, and father GWAS and computes WLM-adjusted effect estimates.
#
##


# Libraries

library(here)
library(conflicted)
library(yaml)
library(glue)
library(janitor)
library(dplyr)

# Solve name space conflicts
conflicts_prefer(dplyr::filter)


# Files

args <- commandArgs(TRUE)

parameters_file <- args[1]

if (!file.exists(parameters_file)) {
  
  stop(glue("Parameters file {parameters_file} not found."))
  
}

look_up_folder <- args[2]

if (!dir.exists(look_up_folder)) {
  
  stop(glue("Regenie folder {pheno_folder} not found."))
  
}

output_file <- args[3]


# Import config file

config <- read_yaml(parameters_file)


# WLM overlap parameters

child_mother_overlap <- 0.2653
child_father_overlap <- 0.197
mother_father_overlap <- 0.0082

# Load files

association_results_merged <- NULL

for (analysis_id in names(config$analyses)) {
  
  analysis <- config$analyses[[analysis_id]]
  
  association_results_analysis <- NULL
  
  for (population in analysis$populations) {
    
    association_results <- read.table(
      file = file.path(look_up_folder, analysis_id, glue("pop_{population}_pheno_{analysis$phenotype}.txt")),
      sep = " ",
      header = T
    ) %>%
      clean_names() %>%
      select(
        chr = chrom,
        pos = genpos,
        id,
        allele0,
        allele1,
        a1freq,
        n,
        beta,
        se,
        log10p
      )
    association_results$analysis <- analysis$name
    association_results$phenotype <- analysis$phenotype
    association_results$population <- population
    association_results$statistics <- "GWAS"
    
    association_results_analysis <- rbind(association_results_analysis, association_results)
    
  }
  
  if ("children" %in% analysis$populations && 
      "mothers" %in% analysis$populations && 
      "fathers" %in% analysis$populations) {
    
    for (variant in unique(association_results_analysis$id)) {
      
      variant_results <- association_results_analysis %>% 
        filter(
          id == variant
        )
      
      if (nrow(variant_results) != 3) {
        
        stop(glue("3 values expected for WLM, {nrow(variant_results)} found for {variant} in {analysis$name}."))
        
      }
      
      beta_child <- variant_results$beta[variant_results$population == "children"]
      beta_mother <- variant_results$beta[variant_results$population == "mothers"]
      beta_father <- variant_results$beta[variant_results$population == "fathers"]
      
      se_child <- variant_results$se[variant_results$population == "children"]
      se_mother <- variant_results$se[variant_results$population == "mothers"]
      se_father <- variant_results$se[variant_results$population == "fathers"]
      
      se_child_2 <- se_child * se_child
      se_mother_2 <- se_mother * se_mother
      se_father_2 <- se_father * se_father
      
      wlm_beta_child <- 2 * beta_child - beta_mother - beta_father
      wlm_se_child <- sqrt(4 * se_child_2 + se_mother_2 + se_father_2 - 4 * child_mother_overlap * sqrt(se_child_2 * se_mother_2) - 4 * child_father_overlap * sqrt(se_child_2 * se_father_2) + 2 * mother_father_overlap * sqrt(se_mother_2 * se_father_2))
      wlm_p_child <- -log10(2 * pnorm(-abs(wlm_beta_child / wlm_se_child)))
      
      wlm_beta_mother <- (3 * beta_mother - 2 * beta_child + beta_father) / 2
      wlm_se_mother <- sqrt(9 * se_mother_2 / 4 + se_child_2 + se_father_2 / 4 - 3 * child_mother_overlap * sqrt(se_child_2 * se_mother_2) - child_father_overlap * sqrt(se_child_2 * se_father_2) + 3 * mother_father_overlap * sqrt(se_mother_2 * se_father_2) / 2)
      wlm_p_mother <- -log10(2 * pnorm(-abs(wlm_beta_mother / wlm_se_mother)))
      
      wlm_beta_father <- (3 * beta_father - 2 * beta_child + beta_mother) / 2
      wlm_se_father <- sqrt(9 * se_father_2 / 4 + se_child_2 + se_mother_2 / 4 - child_mother_overlap * sqrt(se_child_2 * se_mother_2) - 3 * child_father_overlap * sqrt(se_child_2 * se_father_2) + 3 * mother_father_overlap * sqrt(se_mother_2 * se_father_2) / 2)
      wlm_p_father <- -log10(2 * pnorm(-abs(wlm_beta_father / wlm_se_father)))
      
      variant_results_child_wlm <- variant_results %>% filter(population == "children")
      variant_results_child_wlm$beta <- wlm_beta_child
      variant_results_child_wlm$se <- wlm_se_child
      variant_results_child_wlm$log10p <- wlm_p_child
      variant_results_child_wlm$statistics <- "WLM"
      
      variant_results_mother_wlm <- variant_results %>% filter(population == "mothers")
      variant_results_mother_wlm$beta <- wlm_beta_mother
      variant_results_mother_wlm$se <- wlm_se_mother
      variant_results_mother_wlm$log10p <- wlm_p_mother
      variant_results_mother_wlm$statistics <- "WLM"
      
      variant_results_father_wlm <- variant_results %>% filter(population == "fathers")
      variant_results_father_wlm$beta <- wlm_beta_father
      variant_results_father_wlm$se <- wlm_se_father
      variant_results_father_wlm$log10p <- wlm_p_father
      variant_results_father_wlm$statistics <- "WLM"
      
      association_results_analysis <- rbind(association_results_analysis, variant_results_child_wlm)
      association_results_analysis <- rbind(association_results_analysis, variant_results_mother_wlm)
      association_results_analysis <- rbind(association_results_analysis, variant_results_father_wlm)
      
    }
  }
  
  association_results_merged <- rbind(association_results_merged, association_results_analysis)
  
}

# Save to file

write.table(
  x = association_results_merged,
  file = gzfile(output_file),
  col.names = T,
  row.names = F,
  sep = "\t",
  quote = F
)
