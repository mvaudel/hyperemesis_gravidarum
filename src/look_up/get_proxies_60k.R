
##
#
# This script get proxies available in MoBa for the SNPs in the weights tables
#
##

# Libraries

library(conflicted)
library(janitor)
library(yaml)
library(glue)
library(tidyverse)

conflicts_prefer(dplyr::filter)

# File import

config = read_yaml("src/look_up/look_up.yaml")

look_up_data = read.table(
  file = "src/look_up/resources/look_up_merged_23-07-21.gz",
  header = T,
  sep = "\t"
) %>% 
  filter(
    id %in% config$variants
  )

top_hits_table <- read.table(
  file = "consequences/utils/top_hits",
  header = T,
  sep = " "
) %>% 
  clean_names()

meta_results <- read.table(
  file = "consequences/resources/meta",
  header = T,
  sep = "\t"
) %>% 
  clean_names()

# MoBa variants

moba_variants_table <- read.table(
  file = "/home/marc/Projects/MoBa/gwas/MoBaPsychGen_v1-ec-eur-batch-basic-qc.snpstats.gz",
  header = T,
  sep = "\t"
) %>% 
  clean_names() %>% 
  mutate(
    id = paste0(bp, ":", a1, ":", a2),
    id_swap = paste0(bp, ":", a2, ":", a1)
  )


# Check if the variants are in MoBa, if not, find proxies

meta_results_top_hits <- meta_results %>% 
  filter(
    snp %in% top_hits_table$snp
  )

meta_results_top_hits$proxy = NA
meta_results_top_hits$proxy_a0 = NA
meta_results_top_hits$proxy_a1 = NA
meta_results_top_hits$proxy_r2 = NA

missing <- meta_results_top_hits %>% 
  filter(
    !snp %in% look_up_data$id
  ) %>% 
  select(
    chr, pos, snp, a0, a1
  )

for (i in 1:nrow(meta_results_top_hits)) {
  
  id <- meta_results_top_hits$snp[i]
  
  if (id %in% missing$snp) {
    
    chr <- meta_results_top_hits$chr[i]
    id1 <- paste0(meta_results_top_hits$pos[i], ":", meta_results_top_hits$a0[i], ":", meta_results_top_hits$a1[i])
    id2 <- paste0(meta_results_top_hits$pos[i], ":", meta_results_top_hits$a1[i], ":", meta_results_top_hits$a0[i])
    
    print(glue("{Sys.time()}    Loading LD annotation for chromosome {chr}"))
    
    proxies_table <- read.table(
      file = glue("/home/marc/Projects/top_ld/EUR_chr{chr}_no_filter_0.2_1000000_LD.csv.gz"),
      header = T,
      sep = ","
    ) %>% 
      clean_names()
    
    # Find proxies
    
    proxies11 <- proxies_table %>% 
      filter(
        uniq_id_1 == id1
      ) %>% 
      select(
        proxy_id = uniq_id_2,
        r2,
        x_corr
      ) %>% 
      mutate(
        swap = ifelse(x_corr == "+", 0, 1)
      )
    
    proxies12 <- proxies_table %>% 
      filter(
        uniq_id_2 == id1
      ) %>% 
      select(
        proxy_id = uniq_id_1,
        r2,
        x_corr
      ) %>% 
      mutate(
        swap = ifelse(x_corr == "+", 0, 1)
      )
    
    proxies21 <- proxies_table %>% 
      filter(
        uniq_id_1 == id2
      ) %>% 
      select(
        proxy_id = uniq_id_2,
        r2,
        x_corr
      ) %>% 
      mutate(
        swap = ifelse(x_corr == "-", 0, 1)
      )
    
    proxies22 <- proxies_table %>% 
      filter(
        uniq_id_2 == id2
      ) %>% 
      select(
        proxy_id = uniq_id_1,
        r2,
        x_corr
      ) %>% 
      mutate(
        swap = ifelse(x_corr == "-", 0, 1)
      )
    
    proxies <- rbind(proxies11, proxies12, proxies21, proxies22)
    
    if (nrow(proxies) > 0) {
      
      proxies_matched <- proxies %>% 
        filter(
          proxy_id %in% moba_variants_table$id | proxy_id %in% moba_variants_table$id_swap
        ) %>% 
        arrange(
          desc(r2)
        )
      
      if (nrow(proxies_matched) > 0) {
        
        proxy_id <- proxies_matched$proxy_id[1]
        proxy_id_split <- strsplit(proxy_id, ":")[[1]]
        
        j <- which(moba_variants_table$id == proxy_id | moba_variants_table$id_swap == proxy_id)
        
        proxy_moba_id <- moba_variants_table$snp[j]
        
        meta_results_top_hits$proxy[i] <- proxy_moba_id
        meta_results_top_hits$proxy_r2[i] <- proxies_matched$r2[1]
        
        if (proxies_matched$swap[1] == 0) {
          
          meta_results_top_hits$proxy_a0[i] <- proxy_id_split[2]
          meta_results_top_hits$proxy_a1[i] <- proxy_id_split[3]
          
        } else if (proxies_matched$swap[1] == 1) {
          
          meta_results_top_hits$proxy_a0[i] <- proxy_id_split[3]
          meta_results_top_hits$proxy_a1[i] <- proxy_id_split[2]
          
        } else {
          
          stop("Invalid swap value")
          
        }
      }
    }
  }
}

write.table(
  x = meta_results_top_hits,
  file = "consequences/utils/top_hits_proxies",
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\t"
)

ids_alleles <- meta_results_top_hits %>% 
  mutate(
    id = ifelse(is.na(proxy), snp, proxy),
    allele = ifelse(is.na(proxy), a1, proxy_a1)
  ) %>% 
  select(
    id, allele
  )

write.table(
  x = ids_alleles,
  file = "consequences/utils/top_hits_proxies_ids_allele",
  col.names = F,
  row.names = F,
  quote = F,
  sep = " "
)

