
check_preds_full_ranks <- function(preds_df) {
  svd_out <- preds_df %>%
    tibble::column_to_rownames(var="subject_id") %>%
    as.matrix() %>%
    svd()
  stopifnot(!any(dplyr::near(svd_out$d, 0)))
}

get_r2 <- function(target, pred) {
  rss <- sum((target - pred)**2)
  tss <- sum((target - mean(target))**2)
  r2 <- 1 - (rss / tss)
  return(r2)
}

get_mse <- function(target, pred) {
  mse <- mean((target - pred)**2)
  return(mse)
}

get_spearman <- function(target, pred) {
  spearman <- cor(x=target, y=pred, method="spearman")
  return(spearman)
}

get_spearman_pval <- function(target, pred, nperm=1000) {
  true_cor <- cor(x=target, y=pred, method="spearman")
  permut_cor <- 
    purrr::map_dbl(1:nperm, ~ cor(x=target, y=sample(pred), method="spearman"))
  pval <- mean(permut_cor > true_cor)
  return(pval)
}

get_oob_perf <- function(model_df) {
  stopifnot(all(c("target", "subject_id", "baseline") %in% colnames(model_df)))
  
  model <- ranger::ranger(formula = target ~ ., 
                          data = (model_df %>% dplyr::select(-subject_id)), 
                          num.trees=500, importance="permutation")
  tibble::tibble(
    mode = "oob",
    mse = get_mse(model_df$target, model$predictions),
    r2 = get_r2(model_df$target, model$predictions),
    srho = get_spearman(model_df$target, model$predictions)
  )
}

get_loocv_perf <- function(model_df) {
  stopifnot(all(c("target", "subject_id", "baseline") %in% colnames(model_df)))
  
  loocv_out <- purrr::map(sample(1:nrow(model_df)), function(idx) {
    train_df <- model_df[-idx, ] %>% dplyr::select(-subject_id)
    test_df <- model_df[idx, ] %>% dplyr::select(-subject_id)
    model <- ranger::ranger(formula = target ~ ., 
                            data = train_df, 
                            num.trees = 500)
    predictions <- predict(model, test_df)$predictions
    tibble::tibble(idx=idx, target = test_df$target, prediction=predictions)
  }) %>%
    dplyr::bind_rows()
  
  tibble::tibble(
    mode = "loocv",
    mse = get_mse(loocv_out$target, loocv_out$prediction),
    r2 = get_r2(loocv_out$target, loocv_out$prediction),
    srho = get_spearman(loocv_out$target, loocv_out$prediction)
  )
}

get_cross_cohort_perf_combinations <- function(model_df, meta_data) {
  stopifnot(all(c("target", "subject_id", "baseline") %in% colnames(model_df)))
  
  model_dataset_df <- model_df %>%
    dplyr::left_join(meta_data %>% 
                       dplyr::select(subject_id, dataset) %>%
                       dplyr::distinct(),
                     by="subject_id")
  
  all_datasets <- unique(model_dataset_df$dataset)
  
  cross_dataset_out <- purrr::map(all_datasets, function(test_dataset) {
    # test_dataset <- all_datasets[1]
    train_dataset <- all_datasets[!all_datasets %in% test_dataset]
    
    train_df <- model_dataset_df %>%
      dplyr::filter(dataset %in% train_dataset) %>%
      dplyr::select(-c(dataset, subject_id))
    
    test_df <- model_dataset_df %>%
      dplyr::filter(dataset %in% test_dataset) %>%
      dplyr::select(-c(dataset, subject_id))
    
    model <- ranger::ranger(formula = target ~ ., 
                            data = train_df, 
                            num.trees = 500)
    predictions <- predict(model, test_df)$predictions
    tibble::tibble(trainset = paste0(train_dataset, collapse="__"),
                   testset = test_dataset,
                   target = test_df$target, 
                   prediction = predictions,
                   testset_mean = mean(test_df$target))
  }) %>%
    dplyr::bind_rows()
  
  cross_dataset_out %>%
    dplyr::group_by(trainset, testset) %>%
    dplyr::summarise(
      mse = get_mse(target, prediction),
      r2 = get_r2(target, prediction),
      srho = get_spearman(target, prediction),
      mse_same_cohort = get_mse(target, testset_mean),
      r2_same_cohort = get_r2(target, testset_mean),
      .groups = "drop"
    )
}

get_cross_cohort_perf_combinations_repeated <- function(model_df, meta_data, n_iter=25) {
  # TODO
  df_long <- purrr::map(1:n_iter, function(iter) {
    get_cross_cohort_perf_combinations(model_df=model_df, meta_data=meta_data) %>%
      dplyr::mutate(iteration=iter)
  }) %>%
    dplyr::bind_rows()
  df_out <- df_long %>%
    dplyr::group_by(trainset, testset) %>%
    dplyr::summarise(srho_mean = mean(test_srho),
                     srho_sd = sd(test_srho),
                     srho_baseline = mean(srho_baseline),
                     train_n = mean(train_n),
                     test_n = mean(test_n), 
                     .groups = "drop")
  return(df_out)
}


get_cross_cohort_perf_single <- function(model_df, meta_data) {
  stopifnot(all(c("target", "subject_id", "baseline") %in% colnames(model_df)))
  
  model_dataset_df <- model_df %>%
    dplyr::left_join(meta_data %>% 
                       dplyr::select(subject_id, dataset) %>%
                       dplyr::distinct(),
                     by="subject_id")
  
  all_datasets <- unique(model_dataset_df$dataset)
  
  n_per_dataset <- model_dataset_df %>%
    dplyr::count(dataset) %>%
    dplyr::pull(n, dataset)
  
  cross_dataset_out <- purrr::map(all_datasets, function(test_dataset) {
    
    purrr::map( all_datasets[!all_datasets %in% test_dataset], function(train_dataset) {
      
      train_df <- model_dataset_df %>%
        dplyr::filter(dataset %in% train_dataset) %>%
        dplyr::select(-c(dataset, subject_id))
      
      test_df <- model_dataset_df %>%
        dplyr::filter(dataset %in% test_dataset) %>%
        dplyr::select(-c(dataset, subject_id))
      
      model <- ranger::ranger(formula = target ~ ., 
                              data = train_df, 
                              num.trees = 500)
      predictions <- predict(model, test_df)$predictions
      tibble::tibble(trainset = train_dataset,
                     testset = test_dataset,
                     test_baseline = test_df$baseline,
                     target = test_df$target, 
                     prediction = predictions,
                     testset_mean = mean(test_df$target))
    }) %>%
      dplyr::bind_rows()
  }) %>%
    dplyr::bind_rows()
  
  cross_dataset_out %>%
    dplyr::group_by(trainset, testset) %>%
    dplyr::summarise(
      mse = get_mse(target, prediction),
      r2 = get_r2(target, prediction),
      srho = get_spearman(target, prediction),
      srho_baseline = get_spearman(target, test_baseline),
      mse_tmean = get_mse(target, testset_mean),
      .groups = "drop"
    ) %>%
    dplyr::mutate(train_n = n_per_dataset[trainset]) %>%
    dplyr::mutate(test_n = n_per_dataset[testset]) %>%
    dplyr::mutate(trainset = str_remove_all(trainset, "_dataset")) %>%
    dplyr::mutate(testset = str_remove_all(testset, "_dataset"))
}

get_cross_cohort_perf_single_repeated <- function(model_df, meta_data, n_iter=25) {
  df_long <- purrr::map(1:n_iter, function(iter) {
    get_cross_cohort_perf_single(model_df=model_df, meta_data=meta_data) %>%
      dplyr::mutate(iteration=iter)
  }) %>%
    dplyr::bind_rows()
  df_out <- df_long %>%
    dplyr::group_by(trainset, testset) %>%
    dplyr::summarise(srho_mean = mean(srho),
                     srho_sd = sd(srho),
                     srho_baseline = mean(srho_baseline),
                     train_n = mean(train_n),
                     test_n = mean(test_n), 
                     .groups = "drop")
  return(df_out)
  }

get_metadata_covariates <- function(meta_data) {
  specimen_list <- get_specimen_per_day(meta_data)
  
  meta_data_covariates <- meta_data %>%
    dplyr::select(subject_id, infancy_vac, biological_sex, ethnicity, age_at_boost) %>%
    dplyr::distinct() %>%
    dplyr::mutate(infancy_vac_value=1, biological_sex_value=1, ethnicity_value=1) %>%
    tidyr::pivot_wider(names_from=infancy_vac, values_from=infancy_vac_value, 
                       names_prefix="infancy_vac_") %>%
    tidyr::pivot_wider(names_from=biological_sex, values_from=biological_sex_value,
                       names_prefix="biological_sex_") %>%
    tidyr::pivot_wider(names_from=ethnicity, values_from=ethnicity_value, 
                       names_prefix="ethnicity_") %>%
    dplyr::mutate(across(everything(), .fns = ~replace_na(.,0))) %>%
    dplyr::select(-c(infancy_vac_wP, biological_sex_Female, ethnicity_Unknown))
  
  colnames(meta_data_covariates) <- str_replace_all(colnames(meta_data_covariates), " ", "_")
  
  check_preds_full_ranks(meta_data_covariates)
  
  return(meta_data_covariates)
}
