
library(readr)
library(dplyr)
library(stringr)
library(lubridate)

read_gene_meta <- function(input_dir) {
  readr::read_delim(file.path(input_dir, "meta_data", "gene.csv"), 
                    delim=";", show_col_types = FALSE) %>%
    return()
}

read_protein_meta <- function(input_dir) {
  readr::read_delim(file.path(input_dir, "meta_data", "protein.csv"), 
                    delim=";", show_col_types = FALSE) %>%
    return()
}


read_celltype_meta <- function(input_dir) {
  readr::read_delim(file.path(input_dir, "meta_data", "celltype.csv"), 
                    delim=";", show_col_types = FALSE) %>%
    dplyr::mutate(cell_type_name = purrr::map_chr(.data$cell_type_name,
                                                  ~ str_replace_all(.x, " ", "_"))) %>%
    return()
}

read_harmonized_meta_data <- function(input_dir) {
  meta_data <- bind_rows(
    readr::read_tsv(file.path(input_dir, 
                              "harmonized_and_processed_data",
                              "master_harmonized_data_TSV",
                              "training_subject_specimen.tsv"), 
                    show_col_types = FALSE) %>% 
      dplyr::mutate(partition="train"),
    readr::read_tsv(file.path(input_dir, 
                              "harmonized_and_processed_data",
                              "master_harmonized_data_TSV",
                              "challenge_subject_specimen.tsv"), 
                    show_col_types = FALSE) %>% 
      dplyr::mutate(partition="challenge")
  ) %>%
    mutate(age_at_boost = lubridate::interval(ymd(year_of_birth), ymd(date_of_boost)) / years(1)) %>%
    dplyr::mutate(diff_relative_to_boost = (planned_day_relative_to_boost-actual_day_relative_to_boost))
  
  stopifnot(all(meta_data$planned_day_relative_to_boost == meta_data$timepoint))
  stopifnot(all(meta_data$specimen_type == "Blood"))
  
  meta_data <- meta_data %>% 
    dplyr::select(-timepoint) %>% # same as planned_day_relative_to_boost
    dplyr::select(-specimen_type) # all blood
  
  return(meta_data)
}

filter_harmonized_meta_data <- function(meta_data, experimental_data) {
  # 1. Only keep subjects where we have some baseline measurement (planned day 0) before the booster administration
  subject_ids_to_remove <- meta_data %>%
    dplyr::filter(planned_day_relative_to_boost==0) %>%
    dplyr::filter(actual_day_relative_to_boost > 0) %>%
    dplyr::pull(subject_id)
  meta_data <- meta_data %>% dplyr::filter(!(subject_id %in% subject_ids_to_remove))
  
  # 2. Only keep specimen for which we have a least one assay
  assays_per_specimen <- purrr::imap(experimental_data, ~ .x %>% 
                                       dplyr::select(specimen_id) %>%
                                       dplyr::distinct() %>%
                                       dplyr::mutate(assay=.y)) %>%
    dplyr::bind_rows()
  meta_data <- meta_data %>% dplyr::filter(specimen_id %in% assays_per_specimen$specimen_id)
  
  # 3. TODO 
  
  return(meta_data)
}


read_harmonized_experimental_data <- function(input_dir) {
  # input_dir <- "/Users/pschafer/Projects/cmi_pb_third_challenge/data"
  experimental_data <- list.files(path = file.path(input_dir, 
                                                   "harmonized_and_processed_data",
                                                   "master_harmonized_data_TSV"), 
                          pattern = "^training_.*_long.tsv$", full.names = TRUE) %>%
    purrr::set_names() %>%
    purrr::map(~ dplyr::bind_rows(readr::read_tsv(.x, show_col_types = FALSE),
                                  readr::read_tsv(str_replace(.x, "training", "challenge"), show_col_types = FALSE)))
  
  names(experimental_data) <- stringr::str_split(names(experimental_data), "/") %>%
    purrr::map_chr(~ .x[length(.x)]) %>%
    stringr::str_remove("training_") %>%
    stringr::str_remove("_long.tsv")
  
  # TODO: Replace this by my defaults file or not?
  feature_name_per_modality <- list(
    "pbmc_cell_frequency" = "cell_type_name",
    "pbmc_gene_expression" = "versioned_ensembl_gene_id",
    "plasma_antibody_levels" = "isotype_antigen",
    "plasma_cytokine_concentrations_by_legendplex" = "protein_id",
    "plasma_cytokine_concentrations_by_olink" = "protein_id",
    "t_cell_activation" = "stimulation",
    "t_cell_polarization" = "protein_id"
  )
  
  # make sure the feature names have no spaces (such that the wide format works better)
  experimental_data <- purrr::imap(experimental_data, function(df, modality) {
    # modality <- "pbmc_cell_frequency"; df <- experimental_data$pbmc_cell_frequency
    df %>%
      dplyr::mutate(!!feature_name_per_modality[[modality]] := str_replace_all(df[[feature_name_per_modality[[modality]]]], " ", "_"))
  })
  
  return(experimental_data)
}


read_raw_experimental_data <- function(input_dir) {
  # input_dir <- "/Users/pschafer/Projects/cmi_pb_third_challenge/data"
  
  all_exp_files <- list.files(path = file.path(input_dir, "raw_datasets"), 
                              recursive=TRUE, full.names=TRUE) %>%
    magrittr::extract(!grepl("specimen|subject", .))
  
  all_exp_data <- list(
    "pbmc_cell_frequency",
    "pbmc_gene_expression",
    "plasma_ab_titer",
    "plasma_cytokine_concentration_by_legendplex",
    "plasma_cytokine_concentration_by_olink",
    "t_cell_activation",
    "t_cell_polarization"
  ) %>%
    purrr::set_names() %>%
    purrr::map(., function(modality) {
      # modality <- "pbmc_cell_frequency"
      purrr::map(all_exp_files[grepl(modality, all_exp_files)],
                 ~ readr::read_tsv(.x, show_col_types = FALSE)) %>%
        dplyr::bind_rows()
    })
  
  # recreate the isotype_antigen column from the harmonized data
  all_exp_data$plasma_ab_titer <- all_exp_data$plasma_ab_titer %>%
    dplyr::mutate(isotype_antigen = paste0(isotype, "_", antigen))
  
  # TODO: Replace this by my defaults file or not?
  feature_name_per_modality <- list(
    "pbmc_cell_frequency" = "cell_type_name",
    "pbmc_gene_expression" = "versioned_ensembl_gene_id",
    "plasma_ab_titer" = "isotype_antigen",
    "plasma_cytokine_concentration_by_legendplex" = "protein_id",
    "plasma_cytokine_concentration_by_olink" = "protein_id",
    "t_cell_activation" = "stimulation",
    "t_cell_polarization" = "protein_id"
  )
  
  # make sure the feature names have no spaces (such that the wide format works better)
  all_exp_data <- purrr::imap(all_exp_data, function(df, modality) {
    # modality <- "pbmc_cell_frequency"; df <- all_exp_data$pbmc_cell_frequency
    df %>%
      dplyr::mutate(!!feature_name_per_modality[[modality]] := str_replace_all(df[[feature_name_per_modality[[modality]]]], " ", "_"))
  })
  
  return(all_exp_data)
}


filter_harmonized_experimental_data <- function(meta_data, experimental_data, verbose=TRUE) {
  # 1. Only keep specimen that are documented in our meta_data
  experimental_data <- purrr::imap(
    experimental_data, 
    ~ {
      initial_count <- nrow(.x)
      filtered_data <- .x %>% dplyr::filter(specimen_id %in% meta_data$specimen_id)
      final_count <- nrow(filtered_data)
      removed_count <- initial_count - final_count
      if (verbose) message("Removed ", removed_count, " specimens from ", .y, "because missing in meta data")
      return(filtered_data)
    }
  )
  
  # 2. TODO
  
  return(experimental_data)
}

