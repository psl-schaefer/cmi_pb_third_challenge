
library(readr)
library(dplyr)
library(stringr)
library(lubridate)

# this function should ensure that I can use feature names in formulas
# while keeping feature names unique of course
clean_names <- function(orig_features) {
  new_features <- orig_features %>%
    str_replace_all(., "\\+", "_P_") %>%
    str_replace_all(., "\\-", "_N_") %>%
    str_replace_all(., "\\.", "_") %>%
    str_replace_all(., "\\:", "_") %>%
    str_replace_all(., " ", "_") %>%
    str_replace_all(., "/", "_")
  stopifnot(length(unique(new_features)) == length(unique(orig_features)))
  return(new_features)
}

read_gene_meta <- function(input_dir) {
  readr::read_delim(file.path(input_dir, "meta_data", "gene.csv"), 
                    delim=";", show_col_types = FALSE) %>%
    dplyr::mutate(gene_symbol_clean = clean_names(gene_symbol)) %>%
    dplyr::mutate(versioned_ensembl_gene_id_clean = clean_names(versioned_ensembl_gene_id)) %>%
    return()
}

read_gene_meta_plus <- function(input_dir) {
  readr::read_csv(file.path(input_dir, "meta_data", "gene_plus.csv"), 
                  show_col_types = FALSE) %>%
    return()
}

cytokine_uniprot_mapping <- data.frame(
  cytokine = c("CCL8", "IL33", "CXCL12", "OLR1", "IL27", "IL2", "CXCL9", "TGFA", "IL1B", "IL6", "IL4", "TNFSF12",
               "TSLP", "CCL11", "HGF", "FLT3LG", "IL17F", "IL7", "IL13", "IL18", "CCL13", "TNFSF10", "CXCL10",
               "IFNG", "IL10", "CCL19", "TNF", "IL15", "CCL3", "CXCL8", "MMP12", "CSF2", "CSF3", "VEGFA", "IL17C",
               "EGF", "CCL2", "IL17A", "OSM", "CSF1", "CCL4", "CXCL11", "LTA", "CCL7", "MMP1", "IL5"),
  protein_id = c("P80075", "O95760", "P48061", "P78380", "Q8NEV9_Q14213", "P60568", "Q07325", "P01135", "P01584",
                 "P05231", "P05112", "O43508", "Q969D9", "P51671", "P14210", "P49771", "Q96PD4", "P13232", "P35225",
                 "Q14116", "Q99616", "P50591", "P02778", "P01579", "P22301", "Q99731", "P01375", "P40933", "P10147",
                 "P10145", "P39900", "P04141", "P09919", "P15692", "Q9P0M4", "P01133", "P13500", "Q16552", "P13725",
                 "P09603", "P13236", "O14625", "P01374", "P80098", "P03956", "P05113")
)

read_protein_meta <- function(input_dir) {
  readr::read_delim(file.path(input_dir, "meta_data", "protein.csv"), 
                    delim=";", show_col_types = FALSE) %>%
    dplyr::left_join(cytokine_uniprot_mapping, by=c("uniprot_id" = "protein_id")) %>% 
    dplyr::mutate(uniprot_id_clean = clean_names(uniprot_id)) %>%
    dplyr::mutate(cytokine_clean = clean_names(cytokine)) %>%
    return()
}

read_celltype_meta <- function(input_dir) {
  readr::read_delim(file.path(input_dir, "meta_data", "celltype.csv"), 
                    delim=";", show_col_types = FALSE) %>%
    dplyr::mutate(cell_type_name_clean = clean_names(cell_type_name)) %>%
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

filter_meta_data <- function(meta_data, experimental_data) {
  # 1. Only keep subjects where we have at least one baseline measurement (actual day <= 0) before the booster administration
  subject_ids_to_remove <- meta_data %>%
    dplyr::group_by(subject_id) %>%
    dplyr::summarise(min_actual_day_relative_to_boost = min(actual_day_relative_to_boost)) %>%
    dplyr::filter(min_actual_day_relative_to_boost > 0) %>% 
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

experimental_data_settings <- list(
  pbmc_cell_frequency = list(
    feature_col = "cell_type_name",
    value_col = "percent_live_cell",
    feature_subset = c("Monocytes", # no parent, three children: Classical_Monocytes, Intermediate_Monocytes, Non-Classical_Monocytes
                       "Classical_Monocytes", # only parent: Monocytes; no children
                       "Intermediate_Monocytes", # only parent: Monocytes; no children
                       "Non_N_Classical_Monocytes", # only parent: Monocytes; no children
                       "CD8Tcells", # four children: NaiveCD8, TcmCD8, TemCD8, TemraCD8
                       "CD4Tcells", # four children: NaiveCD4, TcmCD4, TemCD4, TemraCD4
                       "Bcells", # only has one child: B cells (CD19+CD3-CD14-CD56-)
                       "NK", # only has one child: CD56high NK cells
                       "pDC", # has no children and no parent
                       "Basophils" # has no children and no parent
    )
  ),
  pbmc_gene_expression = list(
    feature_col = "versioned_ensembl_gene_id",
    value_col = "raw_count" # tpm
  ),
  plasma_ab_titer = list(
    feature_col = "isotype_antigen",
    value_col = "MFI",
    feature_subset = c(
      "IgG_DT", "IgG_FIM2_3", "IgG_OVA", "IgG_TT", "IgG1_DT", "IgG1_FHA", 
      "IgG1_FIM2_3", "IgG1_OVA", "IgG1_PRN", "IgG1_PT", "IgG1_TT", "IgG2_DT", 
      "IgG2_FHA", "IgG2_FIM2_3", "IgG2_OVA", "IgG2_PRN", "IgG2_PT", "IgG2_TT", 
      "IgG3_DT", "IgG3_FHA", "IgG3_FIM2_3", "IgG3_OVA", "IgG3_PRN", "IgG3_PT", 
      "IgG3_TT", "IgG4_DT", "IgG4_FHA", "IgG4_FIM2_3", "IgG4_OVA", "IgG4_PRN", 
      "IgG4_PT", "IgG4_TT", "IgG_FHA", "IgG_PRN", "IgG_PT"
    ),
    unit = "MFI"
  ),
  plasma_cytokine_concentration_by_legendplex = list(
    feature_col = "protein_id",
    value_col = "concentration",
    outliers = c(661, 657, # day 0
                 658, 689, 903, # day 1
                 638, 659, # day 3
                 661) # day 14
  ),
  plasma_cytokine_concentration_by_olink = list(
    feature_col = "protein_id",
    value_col = "concentration",
    feature_subset = c(
      "P01135", "P01374", "P01584", "P03956", "P04141", 
      "P09919", "P13725", "P40933", "P49771", "P78380", 
      "Q16552", "Q8NEV9_Q14213", "Q96PD4", "Q9P0M4", 
      "O14625", "O43508", "O95760", "P01133", "P01375", 
      "P01579", "P02778", "P05112", "P05231", "P09603", 
      "P10145", "P10147", "P13232", "P13236", "P13500", 
      "P14210", "P15692", "P22301", "P35225", "P39900", 
      "P48061", "P50591", "P51671", "P60568", "P80075", 
      "P80098", "Q07325", "Q14116", "Q99616", "Q99731"
      #"Q969D9" # 20% NA fraction
    ),
    unit = "PG/ML"
  ),
  t_cell_activation = list(
    feature_col = "stimulation",
    value_col = "analyte_percentages"
  ),
  t_cell_polarization = list(
    feature_col = "stimulation_protein_id",
    value_col = "analyte_counts",
    # all DMSO features were removed because they are constant
    feature_subset = c(
      "PHA_P01579", "PHA_Q16552", "PHA_P05113", 
      "PT_P01579", "PT_Q16552", "PT_P05113"
    )
  )
)

read_raw_experimental_data <- function(input_dir) {
  # input_dir <- "/Users/pschafer/Projects/cmi_pb_third_challenge/data"
  
  all_exp_files <- list.files(path = file.path(input_dir, "raw_datasets"), 
                              recursive=TRUE, full.names=TRUE) %>%
    magrittr::extract(!grepl("specimen|subject", .))
  
  all_exp_data <- names(experimental_data_settings) %>%
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
  
  # recreate the stimulation_protein_id column from the harmonized data
  all_exp_data$t_cell_polarization <- all_exp_data$t_cell_polarization %>%
    dplyr::mutate(stimulation_protein_id = paste0(stimulation, "_", protein_id))
  
  # make sure the feature names have no spaces and special symbol
  # such that wide format works; and model formulae work
  all_exp_data <- purrr::imap(all_exp_data, function(df, modality) {
    # modality <- "pbmc_cell_frequency"; df <- all_exp_data$pbmc_cell_frequency
    df %>%
      dplyr::mutate(!!experimental_data_settings[[modality]]$feature_col := 
                      clean_names(.data[[experimental_data_settings[[modality]]$feature_col]]))
  })
  
  return(all_exp_data)
}

filter_experimental_data <- function(meta_data, experimental_data, gene_meta, verbose=TRUE) {
  # 1. Only keep specimen that are documented in our meta_data
  experimental_data <- purrr::imap(experimental_data, function(df, modality) {
    initial_count <- nrow(df)
    filtered_data <- df %>% dplyr::filter(specimen_id %in% meta_data$specimen_id)
    final_count <- nrow(filtered_data)
    removed_count <- initial_count - final_count
    if (verbose & (removed_count > 0)) {
      message(modality, " | Removed ", removed_count, " specimens because missing in meta data")
    }
    return(filtered_data)
  })
  cat("\n")
  
  # 2. Only keep features that are in the subset for the modality
  experimental_data <- purrr::imap(experimental_data, function(df, modality) {
    # df <- experimental_data[["plasma_ab_titer"]]; modality <- "plasma_ab_titer"
    # df <- experimental_data[["pbmc_cell_frequency"]]; modality <- "pbmc_cell_frequency"
    if ("feature_subset" %in% names(experimental_data_settings[[modality]])) {
      feature_col <- experimental_data_settings[[modality]]$feature_col
      initial_feature_count <- length(unique(df[[feature_col]]))
      stopifnot(all(
        experimental_data_settings[[modality]][["feature_subset"]] %in% df[[feature_col]]
        ))
      
      df <- df %>%
        dplyr::filter(.data[[feature_col]] %in% experimental_data_settings[[modality]][["feature_subset"]])
      final_feature_count <- length(unique(df[[feature_col]]))
      removed_feature_count <- initial_feature_count - final_feature_count
      
      if (verbose) {
        message(modality, " | Removed ", initial_feature_count - final_feature_count, " features because not in feature subset")
      }
    }
    return(df)
  })
  cat("\n")
  
  # 3. Filter Olink measurements where QC is warn
  initial_feature_count <- nrow(experimental_data$plasma_cytokine_concentration_by_olink)
  
  experimental_data$plasma_cytokine_concentration_by_olink <- 
    experimental_data$plasma_cytokine_concentration_by_olink %>%
    dplyr::filter(quality_control != "Warning")
  
  final_feature_count <- nrow(experimental_data$plasma_cytokine_concentration_by_olink)
  if (verbose) {
    message("plasma_cytokine_concentration_by_olink | Removed ", initial_feature_count - final_feature_count, " features because qc warning")
  }
  cat("\n")
  
  # 4. Only use data with the same unit (MFI for plasma_ab_titer; PG/ML for plasma_cytokine_concentration_by_olink)
  # TODO: This might be a mistake!
  experimental_data <- purrr::imap(experimental_data, function(df, modality) {
    if ("unit" %in% names(experimental_data_settings[[modality]])) {
      initial_feature_count <- nrow(df)
      df <- df %>%
        dplyr::filter(.data[["unit"]] %in% experimental_data_settings[[modality]][["unit"]])
      final_feature_count <- nrow(df)
      if (verbose) {
        message(modality, " | Removed ", initial_feature_count - final_feature_count, " measurements because wrong unit used")
      }
    }
    return(df)
  })
  cat("\n")
  
  # 5. Remove outlier specimen if exist
  experimental_data <- purrr::imap(experimental_data, function(df, modality) {
    if ("outliers" %in% names(experimental_data_settings[[modality]])) {
      df <- df %>%
        dplyr::filter(!specimen_id %in% experimental_data_settings[[modality]][["outliers"]])
      if (verbose) {
        message(modality, " | Removed ", length(experimental_data_settings[[modality]][["outliers"]]), " because specimen is outlier")
      }
    }
    return(df)
  })
  cat("\n")
  
  # 6. Only keep genes in pbmc gex data that have a unique mapping from ensemble id to gene symbol
  ensemble_to_gene <- tibble(versioned_ensembl_gene_id_clean = unique(experimental_data$pbmc_gene_expression$versioned_ensembl_gene_id)) %>%
    dplyr::left_join(gene_meta, by="versioned_ensembl_gene_id_clean") %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::arrange(desc(n), gene_symbol) %>%
    dplyr::filter(n==1) %>%
    dplyr::ungroup()
  
  experimental_data$pbmc_gene_expression <- experimental_data$pbmc_gene_expression %>%
    dplyr::filter(versioned_ensembl_gene_id %in% ensemble_to_gene$versioned_ensembl_gene_id_clean)
  
  return(experimental_data)
}


generate_wide_experimental_data <- function(experimental_data, impute="zero", verbose=TRUE) {
  wide_experimental_data <- purrr::imap(experimental_data, function(df, modality) {
    #df <- experimental_data[[6]]; modality <- names(experimental_data)[6]
    feature_col <- experimental_data_settings[[modality]]$feature_col
    value_col <- experimental_data_settings[[modality]]$value_col
    mtx <- df %>%
      .[, c("specimen_id", feature_col, value_col)] %>%
      tidyr::pivot_wider(names_from=dplyr::all_of(feature_col), 
                         values_from=dplyr::all_of(value_col)) %>%
      tibble::column_to_rownames(var="specimen_id") %>%
      as.matrix()
    
    na_frac <- mean(is.na(mtx))
    if (!is.null(impute) & (na_frac > 0)) {
      if (impute == "zero") {
        mtx[is.na(mtx)] <- 0
        if (verbose) message(modality, " | NA Fraction: ", na_frac, " | Imputed with zeros")
      }
      else if (impute == "median") {
        medians <- matrixStats::colMedians(mtx, na.rm=TRUE)
        for (col in colnames(mtx)) {
          mtx[is.na(mtx[, col]), col] <- medians[col]
        }
      } 
      else if (impute == "mean") {
        means <- matrixStats::colMeans2(mtx, na.rm=TRUE)
        for (col in colnames(mtx)) {
          mtx[is.na(mtx[, col]), col] <- means[col]
        }
      } 
      else {
        stop(paste0(impute, " impute strategey not implemented"))
      }
    }
    return(mtx)
  })
}











### Appendix ### 
read_harmonized_experimental_data_depr <- function(input_dir) {
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
      dplyr::mutate(!!feature_name_per_modality[[modality]] := clean_names(df[[feature_name_per_modality[[modality]]]]))
  })
  
  return(experimental_data)
}

filter_harmonized_experimental_data_depr <- function(meta_data, experimental_data, verbose=TRUE) {
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
