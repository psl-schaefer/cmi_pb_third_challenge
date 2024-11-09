
acceptable_differences <- list(
  day_0 = list(
    "min_diff" = -7,
    "max_diff" = 0
  ),
  day_1 = list(
    "min_diff" = 1,
    "max_diff" = 2
  ),
  day_3 = list(
    "min_diff" = 3,
    "max_diff" = 6
  ),
  day_14 = list(
    "min_diff" = 12,
    "max_diff" = 16
  ),
  day_30 = list(
    "min_diff" = 28,
    "max_diff" = 40
  )
)

get_specimen_per_day <- function(meta_data) {
  specimen_list <- purrr::map(acceptable_differences, function(d_list) {
    meta_data %>%
      dplyr::filter(actual_day_relative_to_boost >= d_list$min_diff) %>%
      dplyr::filter(actual_day_relative_to_boost <= d_list$max_diff) %>%
      dplyr::group_by(subject_id) %>%
      dplyr::slice_min(abs(actual_day_relative_to_boost)) %>%
      dplyr::ungroup() %>%
      dplyr::select(specimen_id, 
                    subject_id, 
                    dataset,
                    actual_day_relative_to_boost, 
                    planned_day_relative_to_boost)
  })
  
  # check that no specimen is considered for two days?
  # not sure this is actually necessary!
  stopifnot(all(purrr::imap(specimen_list, ~ .x %>% dplyr::mutate(day=.y)) %>%
                  dplyr::bind_rows() %>%
                  dplyr::count(specimen_id) %>%
                  dplyr::pull(n)) == 1)
  
  return(specimen_list)
}

# 1.1 Rank the individuals by IgG antibody levels against pertussis toxin (PT) that we detect in plasma 14 days post booster vaccinations.
# 1.2 Rank the individuals by fold change of IgG antibody levels against pertussis toxin (PT) that we detect in plasma 14 days post booster vaccinations compared to titer values at day 0.
# NOTE: Instead of fold change I will use logFC here
generate_targets_task_1 <- function(meta_data, 
                                    experimental_data, 
                                    experimental_data_settings,
                                    specimen_list) {
  
  day_post_booster = "day_14"
  
  targets_task_1 <- dplyr::bind_rows(
    meta_data %>% 
      dplyr::filter(specimen_id %in% specimen_list[["day_0"]][["specimen_id"]]) %>%
      dplyr::mutate(day = "day_0"),
    meta_data %>% 
      dplyr::filter(specimen_id %in% specimen_list[[day_post_booster]][["specimen_id"]]) %>%
      dplyr::mutate(day = day_post_booster)
  ) %>%
    dplyr::left_join(experimental_data$plasma_ab_titer, by="specimen_id") %>%
    dplyr::filter(.data[[experimental_data_settings$plasma_ab_titer$feature_col]] == "IgG_PT") %>%
    dplyr::select(dplyr::all_of(c("subject_id", "day", 
                                  experimental_data_settings$plasma_ab_titer$value_col))) %>%
    tidyr::pivot_wider(names_from="day",
                       values_from=experimental_data_settings$plasma_ab_titer$value_col) %>%
    dplyr::mutate(fc = .data[[day_post_booster]] / .data[["day_0"]]) %>%
    dplyr::mutate(log_fc = log(fc)) %>%
    dplyr::mutate(baseline=day_0, task11=day_14, task12=log_fc) %>%
    tidyr::drop_na()
  
  return(targets_task_1)
}

# 2.1 Rank the individuals by predicted frequency of Monocytes on day 1 post boost after vaccination.
# 2.2 Rank the individuals by fold change of predicted frequency of Monocytes on day 1 post booster vaccination compared to cell frequency values at day 0.
# NOTE: Instead of fold change I will use logFC here
generate_targets_task_2 <- function(meta_data, 
                                    experimental_data, 
                                    experimental_data_settings,
                                    specimen_list) {
  
  day_post_booster = "day_1"
  
  targets_task_2 <- dplyr::bind_rows(
    meta_data %>% 
      dplyr::filter(specimen_id %in% specimen_list[["day_0"]][["specimen_id"]]) %>%
      dplyr::mutate(day = "day_0"),
    meta_data %>% 
      dplyr::filter(specimen_id %in% specimen_list[[day_post_booster]][["specimen_id"]]) %>%
      dplyr::mutate(day = day_post_booster)
  ) %>%
    dplyr::left_join(experimental_data$pbmc_cell_frequency, by="specimen_id") %>%
    dplyr::filter(.data[[experimental_data_settings$pbmc_cell_frequency$feature_col]] == "Monocytes") %>%
    dplyr::select(dplyr::all_of(c("subject_id", "day", 
                                  experimental_data_settings$pbmc_cell_frequency$value_col))) %>%
    tidyr::pivot_wider(names_from="day",
                       values_from=experimental_data_settings$pbmc_cell_frequency$value_col) %>%
    dplyr::mutate(fc = .data[[day_post_booster]] / .data[["day_0"]]) %>%
    dplyr::mutate(log_fc = log(fc)) %>%
    dplyr::mutate(baseline=day_0, task21=day_1, task22=log_fc) %>%
    tidyr::drop_na()
  
  return(targets_task_2)
}

# 3.1 Rank the individuals by predicted gene expression of CCL3 on day 3 post-booster vaccination.
# 3.2 Rank the individuals by fold change of predicted gene expression of CCL3 on day 3 post booster vaccination compared to gene expression values at day 0.
# NOTE: Instead of directly predicting expression, I will predict log expression
# NOTE: Instead of fold change I will use logFC here
generate_targets_task_3 <- function(meta_data, 
                                    experimental_data, 
                                    experimental_data_settings,
                                    specimen_list,
                                    gene_meta) {
  
  experimental_data$pbmc_gene_expression <- experimental_data$pbmc_gene_expression %>%
    dplyr::mutate(specimen_id = as.numeric(specimen_id)) # make sure the join below works
  
  ccl3_ensemble <- experimental_data$pbmc_gene_expression %>%
    dplyr::select("versioned_ensembl_gene_id") %>%
    dplyr::distinct() %>%
    dplyr::left_join(gene_meta, by=c("versioned_ensembl_gene_id"="versioned_ensembl_gene_id_clean")) %>%
    dplyr::filter(gene_symbol=="CCL3") %>%
    dplyr::pull("versioned_ensembl_gene_id")
  stopifnot(length(ccl3_ensemble)==1)
  
  day_post_booster = "day_3"
  
  targets_task_3 <- dplyr::bind_rows(
    meta_data %>% 
      dplyr::filter(specimen_id %in% specimen_list[["day_0"]][["specimen_id"]]) %>%
      dplyr::mutate(day = "day_0"),
    meta_data %>% 
      dplyr::filter(specimen_id %in% specimen_list[[day_post_booster]][["specimen_id"]]) %>%
      dplyr::mutate(day = day_post_booster)
  ) %>%
    dplyr::left_join(experimental_data$pbmc_gene_expression, by="specimen_id") %>%
    dplyr::filter(.data[[experimental_data_settings$pbmc_gene_expression$feature_col]] == ccl3_ensemble) %>%
    dplyr::select(dplyr::all_of(c("subject_id", "day", 
                                  experimental_data_settings$pbmc_gene_expression$value_col))) %>%
    tidyr::pivot_wider(names_from="day",
                       values_from=experimental_data_settings$pbmc_gene_expression$value_col) %>%
    dplyr::mutate(fc = .data[[day_post_booster]] / .data[["day_0"]]) %>%
    dplyr::mutate(log_fc = log(fc)) %>%
    dplyr::mutate(baseline=log(day_0), task31=log(day_3), task32=log_fc) %>%
    tidyr::drop_na()
  
  return(targets_task_3)
}

# 4.1 Rank the individuals based on their Th1/Th2 (IFN-γ/IL-5) polarization ratio on Day 30 post-booster vaccination.
generate_targets_task_4 <- function(meta_data, 
                                    experimental_data, 
                                    experimental_data_settings,
                                    specimen_list,
                                    protein_meta) {
  
  day_post_booster = "day_30"
  
  # https://discuss.cmi-pb.org/t/announcement-bonus-prediction-task-for-the-3rd-public-cmi-pb-challenge/683
  # ratio of IFNG / IL5
  # NOTE: I assume here that "PT" stimulation is correct
  # NOTE: Instead of the ratio I use the log ratio
  targets <- experimental_data$t_cell_polarization %>%
    dplyr::filter(!is.na(protein_id)) %>%
    dplyr::left_join(protein_meta, by=c("protein_id"="uniprot_id")) %>%
    dplyr::filter(stimulation=="PT", cytokine %in% c("IFNG", "IL5")) %>%
    dplyr::select(specimen_id, cytokine, analyte_counts) %>%
    tidyr::pivot_wider(names_from=cytokine, values_from=analyte_counts) %>%
    dplyr::mutate(ratio = IFNG / IL5) %>%
    dplyr::mutate(log_ratio = log(ratio))
  
  day_0 <- meta_data %>% 
    dplyr::filter(specimen_id %in% specimen_list[["day_0"]][["specimen_id"]]) %>%
    dplyr::mutate(day = "day_0") %>%
    dplyr::select(specimen_id, subject_id, day) %>%
    dplyr::left_join(targets, by="specimen_id") %>%
    tidyr::drop_na() %>%
    dplyr::select(subject_id, log_ratio, day)
  
  day_30 <- meta_data %>% 
    dplyr::filter(specimen_id %in% specimen_list[[day_post_booster]][["specimen_id"]]) %>%
    dplyr::mutate(day = day_post_booster) %>%
    dplyr::select(specimen_id, subject_id, day) %>%
    dplyr::left_join(targets, by="specimen_id") %>%
    tidyr::drop_na() %>%
    dplyr::select(subject_id, log_ratio, day)
  
  targets_task_4 <- dplyr::bind_rows(day_0, day_30) %>%
    tidyr::pivot_wider(names_from=day, values_from=log_ratio) %>%
    dplyr::mutate(baseline=day_0, task41=day_30) %>%
    dplyr::filter(!is.na(task41)) %>%
    dplyr::filter(!is.na(baseline)) # there are 3 subjects with missing baseline, which I also remove...
  
  return(targets_task_4)
}

generate_all_targets <- function(meta_data, 
                                 experimental_data, 
                                 experimental_data_settings,
                                 gene_meta,
                                 protein_meta) {
  
  specimen_list <- get_specimen_per_day(meta_data=meta_data)
  
  task_11 <- generate_targets_task_1(meta_data=meta_data, 
                                     experimental_data=experimental_data, 
                                     experimental_data_settings=experimental_data_settings,
                                     specimen_list=specimen_list) %>%
    dplyr::rename(target=task11) %>%
    dplyr::select(subject_id, baseline, target)
    
  task_12 <- generate_targets_task_1(meta_data=meta_data, 
                                     experimental_data=experimental_data, 
                                     experimental_data_settings=experimental_data_settings,
                                     specimen_list=specimen_list) %>%
    dplyr::rename(target=task12) %>% 
    dplyr::select(subject_id, baseline, target)

  task_21 <- generate_targets_task_2(meta_data=meta_data, 
                                     experimental_data=experimental_data, 
                                     experimental_data_settings=experimental_data_settings,
                                     specimen_list=specimen_list) %>%
    dplyr::rename(target=task21) %>%
    dplyr::select(subject_id, baseline, target)
  
  task_22 <- generate_targets_task_2(meta_data=meta_data, 
                                     experimental_data=experimental_data, 
                                     experimental_data_settings=experimental_data_settings,
                                     specimen_list=specimen_list) %>%
    dplyr::rename(target=task22) %>%
    dplyr::select(subject_id, baseline, target)
    
  task_31 <- generate_targets_task_3(meta_data=meta_data, 
                                     experimental_data=experimental_data, 
                                     experimental_data_settings=experimental_data_settings, 
                                     specimen_list=specimen_list,
                                     gene_meta=gene_meta) %>%
    dplyr::rename(target=task31) %>%
    dplyr::select(subject_id, baseline, target)
  
  task_32 <- generate_targets_task_3(meta_data=meta_data, 
                                     experimental_data=experimental_data, 
                                     experimental_data_settings=experimental_data_settings,
                                     specimen_list=specimen_list,
                                     gene_meta=gene_meta) %>%
    dplyr::rename(target=task32) %>%
    dplyr::select(subject_id, baseline, target)
    
  task_41 <- generate_targets_task_4(meta_data=meta_data, 
                                     experimental_data=experimental_data, 
                                     experimental_data_settings=experimental_data_settings,
                                     specimen_list=specimen_list,
                                     protein_meta=protein_meta) %>%
    dplyr::rename(target=task41) %>%
    dplyr::select(subject_id, baseline, target)
  
  return(list(
    task_11 = task_11,
    task_12 = task_12,
    task_21 = task_21,
    task_22 = task_22,
    task_31 = task_31,
    task_32 = task_32,
    task_41 = task_41
  ))
}

plot_targets <- function(target_df) {
  p <- target_df %>%
    tidyr::pivot_longer(cols=starts_with("day_")) %>%
    dplyr::mutate(name = as.numeric(str_remove(name, "day_"))) %>%
    ggplot(aes(x=name, y=value)) +
    geom_point() +
    geom_line(aes(group=subject_id), alpha=0.5) +
    geom_smooth(method="lm", formula=y~x) +
    ggdark::dark_mode(verbose=FALSE) +
    labs(x="Day", y="Target for Task")
  return(p)
}

generate_baseline_task_1 <- function(meta_data, 
                                     experimental_data, 
                                     experimental_data_settings,
                                     specimen_list) {
  
    meta_data %>% 
    dplyr::filter(specimen_id %in% specimen_list[["day_0"]][["specimen_id"]]) %>%
    dplyr::mutate(day = "day_0") %>%
    dplyr::left_join(experimental_data$plasma_ab_titer, by="specimen_id") %>%
    dplyr::filter(.data[[experimental_data_settings$plasma_ab_titer$feature_col]] == "IgG_PT") %>%
    dplyr::mutate(baseline = .data[[experimental_data_settings$plasma_ab_titer$value_col]]) %>%
    dplyr::select(subject_id, baseline) %>%
    return()
}

# 2.1 Rank the individuals by predicted frequency of Monocytes on day 1 post boost after vaccination.
# 2.2 Rank the individuals by fold change of predicted frequency of Monocytes on day 1 post booster vaccination compared to cell frequency values at day 0.
# NOTE: Instead of fold change I will use logFC here
generate_baseline_task_2 <- function(meta_data, 
                                     experimental_data, 
                                     experimental_data_settings,
                                     specimen_list) {
  
    meta_data %>% 
    dplyr::filter(specimen_id %in% specimen_list[["day_0"]][["specimen_id"]]) %>%
    dplyr::mutate(day = "day_0") %>%
    dplyr::left_join(experimental_data$pbmc_cell_frequency, by="specimen_id") %>%
    dplyr::filter(.data[[experimental_data_settings$pbmc_cell_frequency$feature_col]] == "Monocytes") %>%
    dplyr::mutate(baseline = .data[[experimental_data_settings$pbmc_cell_frequency$value_col]]) %>%
    dplyr::select(subject_id, baseline) %>%
    return()
}

# 3.1 Rank the individuals by predicted gene expression of CCL3 on day 3 post-booster vaccination.
# 3.2 Rank the individuals by fold change of predicted gene expression of CCL3 on day 3 post booster vaccination compared to gene expression values at day 0.
# NOTE: Instead of directly predicting expression, I will predict log expression
# NOTE: Instead of fold change I will use logFC here
generate_baseline_task_3 <- function(meta_data, 
                                     experimental_data, 
                                     experimental_data_settings,
                                     specimen_list,
                                     gene_meta) {
  
  experimental_data$pbmc_gene_expression <- experimental_data$pbmc_gene_expression %>%
    dplyr::mutate(specimen_id = as.numeric(specimen_id)) # make sure the join below works
  
  ccl3_ensemble <- experimental_data$pbmc_gene_expression %>%
    dplyr::select("versioned_ensembl_gene_id") %>%
    dplyr::distinct() %>%
    dplyr::left_join(gene_meta, by=c("versioned_ensembl_gene_id"="versioned_ensembl_gene_id_clean")) %>%
    dplyr::filter(gene_symbol=="CCL3") %>%
    dplyr::pull("versioned_ensembl_gene_id")
  stopifnot(length(ccl3_ensemble)==1)
  
  
  meta_data %>% 
    dplyr::filter(specimen_id %in% specimen_list[["day_0"]][["specimen_id"]]) %>%
    dplyr::mutate(day = "day_0") %>%
    dplyr::left_join(experimental_data$pbmc_gene_expression, by="specimen_id") %>%
    dplyr::filter(.data[[experimental_data_settings$pbmc_gene_expression$feature_col]] == ccl3_ensemble) %>%
    dplyr::mutate(baseline = log(.data[[experimental_data_settings$pbmc_gene_expression$value_col]])) %>%
    dplyr::select(subject_id, baseline) %>%
    return()
}

# 4.1 Rank the individuals based on their Th1/Th2 (IFN-γ/IL-5) polarization ratio on Day 30 post-booster vaccination.
generate_baseline_task_4 <- function(meta_data, 
                                     experimental_data, 
                                     experimental_data_settings,
                                     specimen_list,
                                     protein_meta) {

  # https://discuss.cmi-pb.org/t/announcement-bonus-prediction-task-for-the-3rd-public-cmi-pb-challenge/683
  # ratio of IFNG / IL5
  # NOTE: I assume here that "PT" stimulation is correct
  # NOTE: Instead of the ratio I use the log ratio
  targets <- experimental_data$t_cell_polarization %>%
    dplyr::filter(!is.na(protein_id)) %>%
    dplyr::left_join(protein_meta, by=c("protein_id"="uniprot_id")) %>%
    dplyr::filter(stimulation=="PT", cytokine %in% c("IFNG", "IL5")) %>%
    dplyr::select(specimen_id, cytokine, analyte_counts) %>%
    tidyr::pivot_wider(names_from=cytokine, values_from=analyte_counts) %>%
    dplyr::mutate(ratio = IFNG / IL5) %>%
    dplyr::mutate(log_ratio = log(ratio))
  
  meta_data %>% 
    dplyr::filter(specimen_id %in% specimen_list[["day_0"]][["specimen_id"]]) %>%
    dplyr::select(specimen_id, subject_id) %>%
    dplyr::left_join(targets, by="specimen_id") %>%
    dplyr::mutate(baseline=log_ratio) %>%
    dplyr::select(subject_id, baseline) %>%
    tidyr::drop_na() %>%
    return()
}

generate_all_baselines <- function(meta_data, 
                                   experimental_data, 
                                   experimental_data_settings,
                                   gene_meta,
                                   protein_meta) {
  
  specimen_list <- get_specimen_per_day(meta_data=meta_data)
  
  task_11 <- generate_baseline_task_1(meta_data=meta_data, 
                                     experimental_data=experimental_data, 
                                     experimental_data_settings=experimental_data_settings,
                                     specimen_list=specimen_list) %>%
    dplyr::select(subject_id, baseline)
  
  task_12 <- generate_baseline_task_1(meta_data=meta_data, 
                                     experimental_data=experimental_data, 
                                     experimental_data_settings=experimental_data_settings,
                                     specimen_list=specimen_list) %>%
    dplyr::select(subject_id, baseline)
  
  task_21 <- generate_baseline_task_2(meta_data=meta_data, 
                                     experimental_data=experimental_data, 
                                     experimental_data_settings=experimental_data_settings,
                                     specimen_list=specimen_list) %>%
    dplyr::select(subject_id, baseline)
  
  task_22 <- generate_baseline_task_2(meta_data=meta_data, 
                                     experimental_data=experimental_data, 
                                     experimental_data_settings=experimental_data_settings,
                                     specimen_list=specimen_list) %>%
    dplyr::select(subject_id, baseline)
  
  task_31 <- generate_baseline_task_3(meta_data=meta_data, 
                                     experimental_data=experimental_data, 
                                     experimental_data_settings=experimental_data_settings, 
                                     specimen_list=specimen_list,
                                     gene_meta=gene_meta) %>%
    dplyr::select(subject_id, baseline)
  
  task_32 <- generate_baseline_task_3(meta_data=meta_data, 
                                     experimental_data=experimental_data, 
                                     experimental_data_settings=experimental_data_settings,
                                     specimen_list=specimen_list,
                                     gene_meta=gene_meta) %>%
    dplyr::select(subject_id, baseline)
  
  task_41 <- generate_baseline_task_4(meta_data=meta_data, 
                                     experimental_data=experimental_data, 
                                     experimental_data_settings=experimental_data_settings,
                                     specimen_list=specimen_list,
                                     protein_meta=protein_meta) %>%
    dplyr::select(subject_id, baseline)
  
  return(list(
    task_11 = task_11,
    task_12 = task_12,
    task_21 = task_21,
    task_22 = task_22,
    task_31 = task_31,
    task_32 = task_32,
    task_41 = task_41
  ))
}
