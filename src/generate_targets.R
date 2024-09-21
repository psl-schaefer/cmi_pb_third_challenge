
acceptable_differences <- list(
  day_0 = list(
    "d" = 0,
    "min_diff" = -7,
    "max_diff" = 0
  ),
  day_1 = list(
    "d" = 0,
    "min_diff" = 1,
    "max_diff" = 2
  ),
  day_3 = list(
    "d" = 0,
    "min_diff" = 3,
    "max_diff" = 6
  ),
  day_14 = list(
    "d" = 0,
    "min_diff" = 12,
    "max_diff" = 16
  )
)

get_specimen_per_day <- function(meta_data) {
  specimen_list <- purrr::map(acceptable_differences, function(d_list) {
    meta_data %>%
      dplyr::mutate(diff = actual_day_relative_to_boost - d_list$d) %>%
      dplyr::filter(diff >= d_list$min_diff) %>%
      dplyr::filter(diff <= d_list$max_diff) %>%
      dplyr::group_by(subject_id) %>%
      dplyr::slice_min(abs(diff)) %>%
      dplyr::ungroup() %>%
      dplyr::select(specimen_id, 
                    subject_id, 
                    dataset,
                    actual_day_relative_to_boost, 
                    planned_day_relative_to_boost,
                    diff)
  })
  
  # check that no specimen is considered for two days?
  # not sure this is actually necessary!
  stopifnot(all(purrr::imap(specimen_list, ~ .x %>% dplyr::mutate(day=.y)) %>%
                  dplyr::bind_rows() %>%
                  dplyr::count(specimen_id) %>%
                  dplyr::pull(n)) == 1)
  
  return(specimen_list)
}

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
    dplyr::mutate(baseline=day_0, task11=day_14, task12=fc) %>%
    tidyr::drop_na()
  
  return(targets_task_1)
}

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
    dplyr::mutate(baseline=day_0, task21=day_1, task22=fc) %>%
    tidyr::drop_na()
  
  return(targets_task_2)
}

generate_targets_task_3 <- function(meta_data, 
                                    experimental_data, 
                                    experimental_data_settings,
                                    specimen_list,
                                    gene_meta) {
  
  ccl3_ensemble <- experimental_data$pbmc_gene_expression %>%
    dplyr::select("versioned_ensembl_gene_id") %>%
    dplyr::distinct() %>%
    dplyr::left_join(gene_meta, by="versioned_ensembl_gene_id") %>%
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
    dplyr::mutate(baseline=day_0, task31=day_3, task32=fc) %>%
    tidyr::drop_na()
  
  return(targets_task_3)
}

generate_all_targets <- function(meta_data, 
                                 experimental_data, 
                                 experimental_data_settings,
                                 gene_meta) {
  
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
    
  return(list(
    task_11 = task_11,
    task_12 = task_12,
    task_21 = task_21,
    task_22 = task_22,
    task_31 = task_31,
    task_32 = task_32
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