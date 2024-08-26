# look_up

## Look up of top SNPs

Look up of *bona fide* variants in the GWAS files.

### Libraries and settings

``` r
library(conflicted)
library(yaml)
library(janitor)
library(glue)
library(tidyverse)
```

    ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ✔ ggplot2 3.4.4     ✔ purrr   1.0.2
    ✔ tibble  3.2.1     ✔ dplyr   1.1.4
    ✔ tidyr   1.3.1     ✔ stringr 1.5.1
    ✔ readr   2.1.5     ✔ forcats 1.0.0

``` r
library(ggplot2)
library(scico)
library(ggside)
```

    Registered S3 method overwritten by 'ggside':
      method from   
      +.gg   ggplot2

``` r
library(ragg)
library(grid)

# Solve name space conflicts
conflicts_prefer(dplyr::filter)
```

    [conflicted] Will prefer dplyr::filter over any other package.

``` r
# General parameters
theme_set(theme_bw(base_size = 14))
```

### Load data

``` r
# Config

config = read_yaml("look_up.yaml")

# Locus names
variants_details = read.table(
  file = "resources/look_up_snps_24.06.05",
  header = T,
  sep = "\t"
)
loci = variants_details$locus
names(loci) = variants_details$id

# Summary stats
look_up_data = read.table(
  file = "resources/look_up_merged_24-06-05.gz",
  header = T,
  sep = "\t"
) %>% 
  left_join(
    variants_details %>% 
      select(
        id, meta_id, locus
      ), 
      by = "id"
  )
```

### Swap alleles

``` r
ref_analysis = "Nausea vomiting during pregnancy"

for (variant in unique(look_up_data$id)) {
  
  variant_data = look_up_data %>% 
    filter(
      analysis == ref_analysis & id == variant
    ) %>% 
    slice_max(
      log10p, 
      n = 1,
      with_ties = F,
      na_rm = T
      )
  
  if (variant_data$beta < 0) {
    
    look_up_data$beta[look_up_data$id == variant] = -look_up_data$beta[look_up_data$id == variant]
    
  }
  
}
```

### Functions

``` r
get_forest_plot_wlm = function(
  data,
  phenotype_name
) {
  
  variants = sort(unique(data$id))
  names = paste0(loci[variants], " (", variants, ")")
  
  analysis_data = data %>% 
    filter(
      phenotype == phenotype_name & statistics == "WLM" & population %in% c("children", "mothers")
    ) %>% 
    arrange(
      chr, pos
    ) %>% 
    mutate(
      variant_factor = factor(id, levels = variants),
      population_factor = factor(population, levels = rev(c("children", "mothers"))),
      y = as.numeric(variant_factor) + 0.15 * (as.numeric(population_factor) - 2)
    )
  
  levels(analysis_data$population_factor) = rev(c("Child", "Mother"))
  
  plot = ggplot() +
    geom_vline(
      xintercept = 0
    ) +
    geom_segment(
      data = analysis_data,
      mapping = aes(
        x = beta - qnorm(0.975) * se,
        xend = beta + qnorm(0.975) * se,
        y = y,
        yend = y,
        col = population_factor
      ),
      linewidth = 0.5
    ) +
    geom_segment(
      data = analysis_data,
      mapping = aes(
        x = beta - se,
        xend = beta + se,
        y = y,
        yend = y,
        col = population_factor
      ),
      linewidth = 0.8
    ) +
    geom_point(
      data = analysis_data,
      mapping = aes(
        x = beta,
        y = y,
        col = population_factor
      ),
      size = 3
    ) +
    scale_x_continuous(
      name = "beta [95% CI]"
    ) +
    scale_y_continuous(
      breaks = 1:length(variants),
      labels = names,
      limits = c(0.5, length(variants) + 0.5)
    ) +
    scale_color_manual(
      values = rev(c("black", scico(n = 2, palette = "cork", begin = 0.2, end = 0.8)))
    ) + 
    guides(
      colour = guide_legend(
        reverse = T
        )
      ) +
    theme(
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.title = element_blank()
    )
  
  return(plot)
  
}

get_forest_plot_mother = function(
  data,
  phenotype_name
) {
  
  variants = sort(unique(data$id))
  names = paste0(loci[variants], " (", variants, ")")
  
  analysis_data = data %>% 
    filter(
      phenotype == phenotype_name & statistics == "GWAS" & population == "mothers"
    ) %>% 
    arrange(
      chr, pos
    ) %>% 
    mutate(
      variant_factor = factor(id, levels = variants),
      y = as.numeric(variant_factor)
    )
  
  plot = ggplot() +
    geom_vline(
      xintercept = 0
    ) +
    geom_segment(
      data = analysis_data,
      mapping = aes(
        x = beta - qnorm(0.975) * se,
        xend = beta + qnorm(0.975) * se,
        y = y,
        yend = y
      ),
      linewidth = 0.5
    ) +
    geom_segment(
      data = analysis_data,
      mapping = aes(
        x = beta - se,
        xend = beta + se,
        y = y,
        yend = y
      ),
      linewidth = 0.8
    ) +
    geom_point(
      data = analysis_data,
      mapping = aes(
        x = beta,
        y = y
      ),
      size = 3
    ) +
    scale_x_continuous(
      name = "beta [95% CI]"
    ) +
    scale_y_continuous(
      breaks = 1:length(variants),
      labels = names,
      limits = c(0.5, length(variants) + 0.5)
    ) +
    theme(
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.title = element_blank()
    )
  
  return(plot)
  
}
```

### Forest plots

``` r
if (!dir.exists("plots")) {

  dir.create("plots")
  
}

for (phenotype in unique(look_up_data$phenotype)) {
  
  plot = get_forest_plot_wlm(
    data = look_up_data,
    phenotype_name = phenotype
  )
  
  agg_png(
    filename = glue("plots/{make_clean_names(phenotype)}_forest_wlm.png"),
    height = 600,
    width = 900,
    scaling = 1.5
  )
  grid.draw(plot)
  device = dev.off()
  
  plot = get_forest_plot_mother(
    data = look_up_data,
    phenotype_name = phenotype
  )
  
  agg_png(
    filename = glue("plots/{make_clean_names(phenotype)}_forest_mother.png"),
    height = 600,
    width = 900,
    scaling = 1.5
  )
  grid.draw(plot)
  device = dev.off()
  
  grid.draw(plot + ggtitle(phenotype))
  
}
```

![](look_up_files/figure-commonmark/unnamed-chunk-5-1.png)

![](look_up_files/figure-commonmark/unnamed-chunk-5-2.png)

![](look_up_files/figure-commonmark/unnamed-chunk-5-3.png)

![](look_up_files/figure-commonmark/unnamed-chunk-5-4.png)

![](look_up_files/figure-commonmark/unnamed-chunk-5-5.png)

![](look_up_files/figure-commonmark/unnamed-chunk-5-6.png)

![](look_up_files/figure-commonmark/unnamed-chunk-5-7.png)

![](look_up_files/figure-commonmark/unnamed-chunk-5-8.png)

![](look_up_files/figure-commonmark/unnamed-chunk-5-9.png)

![](look_up_files/figure-commonmark/unnamed-chunk-5-10.png)

![](look_up_files/figure-commonmark/unnamed-chunk-5-11.png)

![](look_up_files/figure-commonmark/unnamed-chunk-5-12.png)

![](look_up_files/figure-commonmark/unnamed-chunk-5-13.png)

![](look_up_files/figure-commonmark/unnamed-chunk-5-14.png)

![](look_up_files/figure-commonmark/unnamed-chunk-5-15.png)

![](look_up_files/figure-commonmark/unnamed-chunk-5-16.png)

![](look_up_files/figure-commonmark/unnamed-chunk-5-17.png)

![](look_up_files/figure-commonmark/unnamed-chunk-5-18.png)

![](look_up_files/figure-commonmark/unnamed-chunk-5-19.png)

![](look_up_files/figure-commonmark/unnamed-chunk-5-20.png)

![](look_up_files/figure-commonmark/unnamed-chunk-5-21.png)

![](look_up_files/figure-commonmark/unnamed-chunk-5-22.png)

![](look_up_files/figure-commonmark/unnamed-chunk-5-23.png)

### Time-resolved plots

``` r
nvp_time = data.frame(
  phenotype = c("nausea_vomiting_before_4w", "nausea_vomiting_5w_8w", "nausea_vomiting_9w_12w", "long_term_nausea_vomiting_13w_16w", "long_term_nausea_vomiting_17w_20w", "long_term_nausea_vomiting_21w_24w", "long_term_nausea_vomiting_25w_28w", "long_term_nausea_vomiting_after_29w"),
  week = c("Before 4", "4 to 8", "9 to 12", "13 to 16", "17 to 20", "21 to 24", "25 to 28", "≥29"),
  week_value = c(2, 6, 10.5, 14.5, 18.5, 22.5, 26.5, 30.5)
)

nvp_time_plot = look_up_data %>% 
  filter(
    phenotype %in% nvp_time$phenotype & population %in% c("children", "mothers") & statistics == "WLM"
  ) %>% 
  left_join(
    nvp_time,
    by = "phenotype"
  ) %>% 
  mutate(
    week_factor = factor(week),
    locus = loci[id],
    name = paste0(locus, " (", id, ")"),
    p_category = case_when(
      log10p > 8 ~ "p < 10-8",
      log10p > 6 ~ "p < 10-6",
      log10p > 4 ~ "p < 10-4",
      log10p > 2 ~ "p < 10-2",
      T        ~ "p >= 10-2"
    ),
    p_category_factor = factor(p_category, levels = c("p >= 10-2", "p < 10-2", "p < 10-4", "p < 10-6", "p < 10-8")),
    dx = (as.numeric(factor(population)) - 1.5) / 2
  )

for (variant_id in unique(nvp_time_plot$id)) {
  
  plot_data = nvp_time_plot %>% 
    filter(
      id == variant_id
    )

time_plot = ggplot() +
  geom_hline(
    yintercept = 0,
    col = "black"
  ) +
  geom_segment(
    data = plot_data,
    mapping = aes(
      x = week_value + dx,
      xend = week_value + dx,
      y = beta - qnorm(0.975) * se,
      yend = beta + qnorm(0.975) * se,
      col = population
    ),
    alpha = 0.8,
    size = 1
  ) +
  geom_line(
    data = plot_data,
    mapping = aes(
      x = week_value + dx,
      y = beta,
      col = population
    ),
    size = 1,
    linetype = "dotted"
  ) +
  geom_point(
    data = plot_data,
    mapping = aes(
      x = week_value + dx,
      y = beta,
      col = population,
      size = p_category_factor
    )
  ) +
  scale_x_continuous(
    name = "gestational week",
    breaks = nvp_time$week_value,
    labels = nvp_time$week
  ) +
  scale_y_continuous(
    name = "beta [95% CI]"
  ) +
  scale_color_manual(
      values = c("black", scico(n = 2, palette = "cork", begin = 0.2, end = 0.8)[1])
  ) +
  scale_size_discrete(
    drop = F
  ) +
    guides(
      colour = guide_legend(
        reverse = F
        ),
      shape = guide_legend(
        reverse = T
        )
      ) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    legend.box = "vertical",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
  
  agg_png(
    filename = glue("plots/{make_clean_names(variant_id)}_nvp_time.png"),
    height = 600,
    width = 900,
    scaling = 1.5
  )
  grid.draw(time_plot)
  device = dev.off()
  
  grid.draw(time_plot + ggtitle(plot_data$name[1] ))
  
}
```

    Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ℹ Please use `linewidth` instead.

    Warning: Using size for a discrete variable is not advised.
    Using size for a discrete variable is not advised.

![](look_up_files/figure-commonmark/unnamed-chunk-6-1.png)

    Warning: Using size for a discrete variable is not advised.

![](look_up_files/figure-commonmark/unnamed-chunk-6-2.png)

    Warning: Using size for a discrete variable is not advised.

![](look_up_files/figure-commonmark/unnamed-chunk-6-3.png)

    Warning: Using size for a discrete variable is not advised.

![](look_up_files/figure-commonmark/unnamed-chunk-6-4.png)

    Warning: Using size for a discrete variable is not advised.

![](look_up_files/figure-commonmark/unnamed-chunk-6-5.png)

    Warning: Using size for a discrete variable is not advised.

![](look_up_files/figure-commonmark/unnamed-chunk-6-6.png)

    Warning: Using size for a discrete variable is not advised.

![](look_up_files/figure-commonmark/unnamed-chunk-6-7.png)

    Warning: Using size for a discrete variable is not advised.

![](look_up_files/figure-commonmark/unnamed-chunk-6-8.png)

![](look_up_files/figure-commonmark/unnamed-chunk-6-9.png)

### All variants, mothers only

``` r
variants_to_plot = variants_details$id
loci_to_plot = variants_details$locus


nvp_others_time_plot = look_up_data %>% 
  filter(
    phenotype %in% nvp_time$phenotype & population %in% c("mothers") & statistics == "GWAS" & id %in% variants_to_plot
  ) %>% 
  left_join(
    nvp_time,
    by = "phenotype"
  ) %>% 
  mutate(
    week_factor = factor(week),
    locus = loci[id],
    name = paste0(locus, " (", id, ")"),
    p_category = case_when(
      log10p > 8 ~ "p < 10-8",
      log10p > 6 ~ "p < 10-6",
      log10p > 4 ~ "p < 10-4",
      log10p > 2 ~ "p < 10-2",
      T        ~ "p >= 10-2"
    ),
    p_category_factor = factor(p_category, levels = c("p >= 10-2", "p < 10-2", "p < 10-4", "p < 10-6", "p < 10-8")),
    id_factor = factor(id, levels = variants_to_plot),
    dx = (as.numeric(id_factor) - length(variants_to_plot) / 2 - 0.5) / length(variants_to_plot)
  )

time_plot_mothers = ggplot() +
  geom_hline(
    yintercept = 0,
    col = "black"
  ) +
  geom_segment(
    data = nvp_others_time_plot,
    mapping = aes(
      x = week_value + dx,
      xend = week_value + dx,
      y = beta - qnorm(0.975) * se,
      yend = beta + qnorm(0.975) * se,
      col = id_factor
    ),
    alpha = 0.8,
    size = 1
  ) +
  geom_line(
    data = nvp_others_time_plot,
    mapping = aes(
      x = week_value + dx,
      y = beta,
      col = id_factor
    ),
    size = 1,
    linetype = "dotted"
  ) +
  geom_point(
    data = nvp_others_time_plot,
    mapping = aes(
      x = week_value + dx,
      y = beta,
      col = id_factor,
      size = p_category_factor
    )
  ) +
  scale_x_continuous(
    name = "gestational week",
    breaks = nvp_time$week_value,
    labels = nvp_time$week
  ) +
  scale_y_continuous(
    name = "beta [95% CI]"
  ) +
  scale_color_manual(
    values = scico(
      n = length(variants_to_plot),
      palette = "batlow",
      begin = 0.2,
      end = 0.8
    ),
    labels = loci_to_plot
  ) +
  scale_size_discrete(
    drop = F
  ) +
    guides(
      colour = guide_legend(
        reverse = F
        ),
      shape = guide_legend(
        reverse = T
        )
      ) +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    legend.box = "vertical",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  facet_grid(
    id_factor ~ .
  )
```

    Warning: Using size for a discrete variable is not advised.

``` r
  agg_png(
    filename = glue("plots/nvp_time_mothers.png"),
    height = 8 * 600,
    width = 900,
    scaling = 1.5
  )
  grid.draw(time_plot_mothers)
  device = dev.off()
  
  time_plot_mothers
```

![](look_up_files/figure-commonmark/unnamed-chunk-7-1.png)

### Table with summary stats for Marlena

``` r
table_marlena <- look_up_data %>% 
  select(
    chr, pos, allele0, allele1, a1freq, id, locus, phenotype, population, statistics, n, beta, se, log10p
  )

write.table(
  table_marlena,
  file = "resources/table_marlena_24.06.05",
  col.names = T,
  row.names = F,
  sep = "\t",
  quote = F
)
```

### Heatmap

``` r
stats <- "GWAS"

locus_max <- look_up_data %>% 
  filter(
    phenotype %in% nvp_time$phenotype & population == "mothers" & statistics == stats
  ) %>% 
  group_by(
    id
  ) %>% 
  summarise(
    max_beta = max(beta)
  )

locus_time <- look_up_data %>% 
  filter(
    phenotype %in% nvp_time$phenotype & population == "mothers" & statistics == stats
  ) %>% 
  left_join(
    nvp_time,
    by = "phenotype"
  ) %>% 
  left_join(
    locus_max,
    by = "id"
  ) %>% 
  filter(
    beta > max_beta / 2
  ) %>% 
  group_by(
    id
  ) %>% 
  summarize(
    start = min(week_value),
    end = max(week_value)
  ) %>% 
  arrange(
    start, end
  ) %>% 
  left_join(
    variants_details %>% 
      select(id, locus),
    by = "id"
  )

locus_order <- c("GDF15", "PGR-TRPC6", "IGFBP7", "GFRAL", "SYN3", "SLITRK1", "IGSF11", "TCF7L2-HABP2", "FSHB")

heatmap_data = look_up_data %>% 
  filter(
    phenotype %in% nvp_time$phenotype & population == "mothers" & statistics == stats
  ) %>% 
  left_join(
    nvp_time,
    by = "phenotype"
  ) %>% 
  mutate(
    week_factor = factor(week, levels = nvp_time$week),
    locus_factor = factor(locus, levels = rev(locus_order)),
    p_category = case_when(
      log10p > 6 ~ "p < 10-6",
      log10p > 4 ~ "p < 10-4",
      log10p > 2 ~ "p < 10-2",
      T        ~ "p >= 10-2"
    ),
    p_category_factor = factor(p_category, levels = c("p >= 10-2", "p < 10-2", "p < 10-4", "p < 10-6"))
  )

names(heatmap_data$locus_factor) <- rev(paste0(locus_time$locus, " (", locus_time$id, ")"))

min_fill = min(heatmap_data$beta)

if (min_fill > 0) {
  
  stop("Min beta expected to be negative")
  
}

max_fill = max(heatmap_data$beta)

if (max_fill < 0) {
  
  stop("Max beta expected to be positive")
  
}

fill_scale = 0.4
abs_fill = max(abs(max_fill), abs(min_fill))
begin_fill = 0.5 + fill_scale * 0.5 * min_fill / abs_fill
end_fill = 0.5 + fill_scale * 0.5 * max_fill / abs_fill

heatmap_gwas <- ggplot(
  data = heatmap_data
) +
  theme_bw(
    base_size = 24
  ) +
  geom_tile(
    mapping = aes(
      x = week_factor,
      y = locus_factor,
      fill = beta
    )
  ) +
  geom_text(
    mapping = aes(
      x = week_factor,
      y = locus_factor,
      label = round(beta, 2),
      col = p_category_factor
    ),
      size = 6
  ) +
  scale_x_discrete(
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    expand = c(0, 0)
  ) +
  scale_fill_scico(
    name = "Beta",
    palette = "vik",
    begin = begin_fill,
    end = end_fill,
    guide = guide_colorbar(
    barwidth = 10
    )
  ) +
  scale_color_manual(
    name = "p-value",
    values = c("white", "grey60", "grey40", "grey20"),
        labels = c("", "", "", "", ""),
        guide = guide_legend(
            override.aes = list(
                label = c("p >= 10-2", "p < 10-2", "p < 10-4", "p < 10-6"),
                size = 6
            )
        )
  ) +
  theme(
    axis.title = element_blank(),
    legend.position = "top",
    legend.spacing.x = unit(1.0, "cm")
  )

png(
  file = "plots/heatmap_gwas.png",
  height = 600,
  width = 900
)
grid.draw(heatmap_gwas)
device <- dev.off()
```

#### Bubble plot

``` r
bubble_gwas <- ggplot(
  data = heatmap_data
) +
  theme_bw(
    base_size = 24
  ) +
  geom_point(
    mapping = aes(
      x = week_factor,
      y = locus_factor,
      col = beta,
      size = p_category_factor
    )
  ) +
  scale_x_discrete(
    name = "Gestational week"
  ) +
  scale_size_manual(
    name = NULL,
    values = c(5, 10, 15, 20)
  ) +
  scale_color_scico(
    name = "Beta",
    palette = "vik",
    begin = begin_fill,
    end = end_fill,
    guide = guide_colorbar(
    barwidth = 12
    )
  ) +
  theme(
    axis.title.y = element_blank(),
    legend.position = "top",
    legend.margin = margin(t = 0, b = 0, l = -40, r = 20, unit = 'pt'),
    panel.grid = element_blank()
  )

png(
  file = "plots/bubble_gwas_24.08.26.png",
  height = 600,
  width = 900
)
grid.draw(bubble_gwas)
device <- dev.off()
```
