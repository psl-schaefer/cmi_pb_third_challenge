
library(tidyverse)
source(file.path("src", "read_data.R"))
source(file.path("src", "generate_targets.R"))

gene_meta <- read_gene_meta(file.path("data"))

meta_data <- read_harmonized_meta_data(file.path("data"))
specimen_per_day <- get_specimen_per_day(meta_data=meta_data)

raw_experimental_data <- read_raw_experimental_data(file.path("data"))
raw_experimental_data <- filter_experimental_data(meta_data=meta_data, 
                                                  experimental_data=raw_experimental_data,
                                                  gene_meta=gene_meta)

gex_df <- raw_experimental_data$pbmc_gene_expression %>%
  dplyr::left_join((meta_data %>% dplyr::select(specimen_id, dataset)),
                   by="specimen_id")
rm(raw_experimental_data)

gex_df_overall <- gex_df %>%
  dplyr::group_by(dataset, versioned_ensembl_gene_id) %>%
  dplyr::summarise(mean_tpm = mean(tpm), var_tpm = var(tpm), .groups="drop_last") %>%
  dplyr::mutate(rank = rank(mean_tpm, na.last=TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(disp = var_tpm / mean_tpm) %>%
  dplyr::arrange(versioned_ensembl_gene_id, dataset)

# gex_df_overall %>%
#   dplyr::filter(dataset=="2023_dataset") %>%
#   ggplot(aes(x=log10(mean_tpm), y=log10(var_tpm))) +
#   geom_point() +
#   geom_smooth(method="lm")

gex_df_baseline <- gex_df %>%
  dplyr::filter(specimen_id %in% specimen_per_day$day_0$specimen_id) %>%
  dplyr::group_by(dataset, versioned_ensembl_gene_id) %>%
  dplyr::summarise(mean_tpm = mean(tpm), var_tpm = var(tpm), .groups="drop_last") %>%
  dplyr::mutate(rank = rank(mean_tpm, na.last=TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(disp = var_tpm / mean_tpm) %>%
  dplyr::arrange(versioned_ensembl_gene_id, dataset)

# check correlations
rank_corr_overall <- purrr::map(unique(gex_df_overall$dataset), function(d1) {
  purrr::map(unique(gex_df_overall$dataset), function(d2) {
    tibble(d1=d1, 
           d2=d2, 
           rho = cor(x=gex_df_overall$mean_tpm[gex_df_overall$dataset==d1],
                     y=gex_df_overall$mean_tpm[gex_df_overall$dataset==d2], 
                     method="pearson"),
           srho = cor(x=gex_df_overall$mean_tpm[gex_df_overall$dataset==d1],
                      y=gex_df_overall$mean_tpm[gex_df_overall$dataset==d2], 
                      method="spearman")
    )
  }) %>%
    dplyr::bind_rows()
}) %>%
  dplyr::bind_rows()

rank_corr_baseline <- purrr::map(unique(gex_df_overall$dataset), function(d1) {
  purrr::map(unique(gex_df_overall$dataset), function(d2) {
    tibble(d1=d1, 
           d2=d2, 
           rho = cor(x=gex_df_baseline$mean_tpm[gex_df_baseline$dataset==d1],
                     y=gex_df_baseline$mean_tpm[gex_df_baseline$dataset==d2], 
                     method="pearson"),
           srho = cor(x=gex_df_baseline$mean_tpm[gex_df_baseline$dataset==d1],
                      y=gex_df_baseline$mean_tpm[gex_df_baseline$dataset==d2], 
                      method="spearman")
    )
  }) %>%
    dplyr::bind_rows()
}) %>%
  dplyr::bind_rows()

gex_df_overall_summary <- gex_df_overall %>%
  dplyr::group_by(versioned_ensembl_gene_id) %>%
  dplyr::summarise(mean_disp = mean(disp),
                   median_disp = median(disp),
                   mean_rank = mean(rank), 
                   mean_tpm = mean(mean_tpm),
                   .groups = "drop")

gene_meta_plus <- gene_meta %>%
  dplyr::left_join(gex_df_overall_summary, 
                   by=c("versioned_ensembl_gene_id_clean"="versioned_ensembl_gene_id"))

readr::write_csv(x = gene_meta_plus,
                 file = file.path("data", "meta_data", "gene_plus.csv"))
