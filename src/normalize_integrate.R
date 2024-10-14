library(tidyverse)
library(DESeq2)
library(sva)


normalize_experimental_data <- function(meta_data, 
                                        raw_experimental_data, 
                                        gene_meta) {
  specimen_per_day <- get_specimen_per_day(meta_data=meta_data)
  normalized_data <- list()
  
  # 1) pbmc_cell_frequency: Median Baseline Normalization
  baseline_medians <- raw_experimental_data$pbmc_cell_frequency %>%
    dplyr::inner_join(specimen_per_day$day_0, by="specimen_id") %>%
    dplyr::group_by(dataset, cell_type_name) %>%
    dplyr::summarise(median = median(percent_live_cell, na.rm=TRUE), 
                     .groups="drop") %>%
    dplyr::mutate(median = ifelse(dplyr::near(median, 0), 1, median)) # prevent division by 0?
  
  normalized_data$pbmc_cell_frequency <- raw_experimental_data$pbmc_cell_frequency %>%
    dplyr::left_join((meta_data %>% dplyr::select(specimen_id, dataset)),
                     by="specimen_id") %>%
    dplyr::left_join(baseline_medians, by=c("dataset", "cell_type_name")) %>%
    dplyr::mutate(percent_live_cell = percent_live_cell / median) %>%
    dplyr::select(-c(dataset, median))
  
  # 2) pbmc_gene_expression: VST
  # a) filter for meaningful genes
  feature_col <- experimental_data_settings$pbmc_gene_expression$feature_col
  value_col <- experimental_data_settings$pbmc_gene_expression$value_col
  pbmc_gene_expression_wide <- raw_experimental_data$pbmc_gene_expression %>%
    dplyr::select(dplyr::all_of(c("specimen_id", feature_col, value_col))) %>%
    tidyr::pivot_wider(names_from=dplyr::all_of(feature_col), 
                       values_from=dplyr::all_of(value_col)) %>%
    tibble::column_to_rownames(var="specimen_id") %>%
    as.matrix()
  
  min_counts <- 200
  min_samples <- 20
  
  # only keep genes that have at least 200 counts; 
  # and are non-zero in at least 20 samples
  genes_to_keep <- colnames(pbmc_gene_expression_wide)[
    (colSums(pbmc_gene_expression_wide) >= min_counts) &
      (colSums(pbmc_gene_expression_wide > 0) >= min_samples)
  ]
  pbmc_gene_expression_wide <- pbmc_gene_expression_wide[ , genes_to_keep]
  
  # check that all column (gene) names are unique
  stopifnot(length(colnames(pbmc_gene_expression_wide)) == length(unique(colnames(pbmc_gene_expression_wide))))
  
  # DESeq2 normalization with variance stabilizing transformation (vst)
  # transforms the count data (normalized by division by the size factors or normalization factors), 
  # yielding a matrix of values which are now approximately homoskedastic
  # The transformation also normalizes with respect to library size
  # perform the transformation per cohort
  pbmc_gene_expression_meta <- tibble::column_to_rownames(meta_data, var="specimen_id") %>%
    .[rownames(pbmc_gene_expression_wide), ] 
  pbmc_gene_expression_normalized <- matrix(data=0, 
                                            nrow=nrow(pbmc_gene_expression_wide),
                                            ncol=ncol(pbmc_gene_expression_wide))
  dimnames(pbmc_gene_expression_normalized) <- dimnames(pbmc_gene_expression_wide)
  
  for (year in unique(pbmc_gene_expression_meta$dataset)) {
    # year <- unique(pbmc_gene_expression_meta$dataset)[1]
    specimen_per_year <- rownames(pbmc_gene_expression_meta)[pbmc_gene_expression_meta$dataset==year]
    normalized_per_year <- 
      DESeq2::DESeqDataSetFromMatrix(countData = t(pbmc_gene_expression_wide[specimen_per_year, ]),
                                     colData = tibble::tibble(specimen_id=specimen_per_year),
                                     design = ~ 1)
    normalized_per_year <- DESeq2::estimateSizeFactors(normalized_per_year, quiet=TRUE) 
    normalized_per_year <- DESeq2::estimateDispersions(normalized_per_year, quiet=TRUE)
    normalized_per_year <- DESeq2::varianceStabilizingTransformation(normalized_per_year, blind=TRUE)
    normalized_per_year <- t(assay(normalized_per_year))
    pbmc_gene_expression_normalized[rownames(normalized_per_year), colnames(normalized_per_year)] <-
      normalized_per_year
  }
  # convert back to long format
  normalized_data$pbmc_gene_expression <- pbmc_gene_expression_normalized %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="specimen_id") %>%
    tidyr::pivot_longer(cols=!c(specimen_id), 
                        names_to=feature_col,
                        values_to=value_col)
  
  normalized_data$pbmc_gene_expression_counts <- pbmc_gene_expression_wide %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="specimen_id") %>%
    tidyr::pivot_longer(cols=!c(specimen_id), 
                        names_to=feature_col,
                        values_to=value_col)
    
  
  # 3) plasma_ab_titer: Median Baseline Normalization
  baseline_medians <- raw_experimental_data$plasma_ab_titer %>%
    dplyr::inner_join(specimen_per_day$day_0, by="specimen_id") %>%
    dplyr::group_by(dataset, isotype_antigen) %>%
    dplyr::summarise(median = median(MFI, na.rm=TRUE), .groups="drop") %>%
    dplyr::mutate(median = ifelse(dplyr::near(median, 0), 1, median)) # prevent division by 0?
  
  normalized_data$plasma_ab_titer <- raw_experimental_data$plasma_ab_titer %>%
    dplyr::left_join((meta_data %>% dplyr::select(specimen_id, dataset)),
                     by="specimen_id") %>%
    dplyr::left_join(baseline_medians, by=c("dataset", "isotype_antigen")) %>%
    dplyr::mutate(MFI = MFI / median) %>%
    dplyr::select(-c(dataset, median))
  
  # 4) plasma_cytokine_concentration_by_legendplex: Median Baseline Normalization
  baseline_medians <- raw_experimental_data$plasma_cytokine_concentration_by_legendplex %>%
    dplyr::inner_join(specimen_per_day$day_0, by="specimen_id") %>%
    dplyr::group_by(dataset, protein_id) %>%
    dplyr::summarise(median = median(concentration, na.rm=TRUE), 
                     .groups="drop") %>%
    dplyr::mutate(median = ifelse(dplyr::near(median, 0), 1, median)) # prevent division by 0?
  
  normalized_data$plasma_cytokine_concentration_by_legendplex <- 
    raw_experimental_data$plasma_cytokine_concentration_by_legendplex %>%
    dplyr::left_join((meta_data %>% dplyr::select(specimen_id, dataset)),
                     by="specimen_id") %>%
    dplyr::left_join(baseline_medians, by=c("dataset", "protein_id")) %>%
    dplyr::mutate(concentration = concentration / median) %>%
    dplyr::select(-c(dataset, median))
  
  # 5) plasma_cytokine_concentration_by_olink: No Normalization needed? 
  # TODO: Take into consideration that the CMI-PB consortium suggests median baseline normalization
  normalized_data$plasma_cytokine_concentration_by_olink <- 
    raw_experimental_data$plasma_cytokine_concentration_by_olink
  
  # 6) t_cell_activation: No Normalization needed? (suggested by CMI-PB consortium)
  normalized_data$t_cell_activation <- 
    raw_experimental_data$t_cell_activation
  
  # 7) t_cell_polarization: No Normalization needed? (suggested by CMI-PB consortium)
  normalized_data$t_cell_polarization <- 
    raw_experimental_data$t_cell_polarization
  
  # check that we did not forget any assay
  stopifnot(dplyr::setequal(names(raw_experimental_data), 
                            names(normalized_data)[!names(normalized_data) %in% c("pbmc_gene_expression_counts")]))
  
  return(normalized_data)
}


integrate_experimental_data <- function(meta_data, normalized_experimental_data) {
  integrated_data <- list()
  
  # 1) pbmc_cell_frequency: Nothing
  integrated_data$pbmc_cell_frequency <-
    normalized_experimental_data$pbmc_cell_frequency
  
  # 2) pbmc_gene_expression: ComBat-seq
  feature_col <- experimental_data_settings$pbmc_gene_expression$feature_col
  value_col <- experimental_data_settings$pbmc_gene_expression$value_col
  pbmc_gene_expression_wide <- normalized_experimental_data$pbmc_gene_expression_counts %>%
    dplyr::select(dplyr::all_of(c("specimen_id", feature_col, value_col))) %>%
    tidyr::pivot_wider(names_from=dplyr::all_of(feature_col), 
                       values_from=dplyr::all_of(value_col)) %>%
    tibble::column_to_rownames(var="specimen_id") %>%
    as.matrix()
  
  pbmc_gene_expression_meta <- tibble::column_to_rownames(meta_data, var="specimen_id") %>%
    .[rownames(pbmc_gene_expression_wide), ] 
  
  pbmc_gene_expression_corrected <- sva::ComBat_seq(
    counts = t(pbmc_gene_expression_wide), batch = pbmc_gene_expression_meta$dataset
    )
  pbmc_gene_expression_corrected <- t(pbmc_gene_expression_corrected)
  pbmc_gene_expression_corrected <- matrix(as.integer(pbmc_gene_expression_corrected), 
                                           nrow = nrow(pbmc_gene_expression_corrected), 
                                           ncol = ncol(pbmc_gene_expression_corrected), 
                                           dimnames = dimnames(pbmc_gene_expression_corrected))
  
  # just to make sure, if any specimen have been removed by combat-seq
  pbmc_gene_expression_meta <- tibble::column_to_rownames(meta_data, var="specimen_id") %>%
    .[rownames(pbmc_gene_expression_corrected), ] 
  
  # again vst normalize
  pbmc_gene_expression_normalized <- matrix(data=0, 
                                            nrow=nrow(pbmc_gene_expression_corrected),
                                            ncol=ncol(pbmc_gene_expression_corrected))
  dimnames(pbmc_gene_expression_normalized) <- dimnames(pbmc_gene_expression_corrected)
  
  for (year in unique(pbmc_gene_expression_meta$dataset)) {
    # year <- unique(pbmc_gene_expression_meta$dataset)[1]
    specimen_per_year <- rownames(pbmc_gene_expression_meta)[pbmc_gene_expression_meta$dataset==year]
    normalized_per_year <- 
      DESeq2::DESeqDataSetFromMatrix(countData = t(pbmc_gene_expression_corrected[specimen_per_year, ]),
                                     colData = tibble::tibble(specimen_id=specimen_per_year),
                                     design = ~ 1)
    normalized_per_year <- DESeq2::estimateSizeFactors(normalized_per_year, quiet=TRUE) 
    normalized_per_year <- DESeq2::estimateDispersions(normalized_per_year, quiet=TRUE)
    normalized_per_year <- DESeq2::varianceStabilizingTransformation(normalized_per_year, blind=TRUE)
    normalized_per_year <- t(assay(normalized_per_year))
    pbmc_gene_expression_normalized[rownames(normalized_per_year), colnames(normalized_per_year)] <-
      normalized_per_year
  }
  # convert back to long format
  integrated_data$pbmc_gene_expression <- 
    pbmc_gene_expression_normalized %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="specimen_id") %>%
    tidyr::pivot_longer(cols=!c(specimen_id), 
                        names_to=feature_col,
                        values_to=value_col)
  
  integrated_data$pbmc_gene_expression_counts <- 
    pbmc_gene_expression_corrected %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="specimen_id") %>%
    tidyr::pivot_longer(cols=!c(specimen_id), 
                        names_to=feature_col,
                        values_to=value_col)
  
  # 3) plasma_ab_titer: 
  # TODO: Note sure yet
  integrated_data$plasma_ab_titer <- 
    normalized_experimental_data$plasma_ab_titer
  
  # 4) plasma_cytokine_concentration_by_legendplex: 
  #integrated_data$plasma_cytokine_concentration_by_legendplex <- 
  #  experimental_data$plasma_cytokine_concentration_by_legendplex
  feature_col <- experimental_data_settings$plasma_cytokine_concentration_by_legendplex$feature_col
  value_col <- experimental_data_settings$plasma_cytokine_concentration_by_legendplex$value_col

  plasma_cytokine_concentration_by_legendplex_wide <- 
    normalized_experimental_data$plasma_cytokine_concentration_by_legendplex %>%
    dplyr::select(dplyr::all_of(c("specimen_id", feature_col, value_col))) %>%
    tidyr::pivot_wider(names_from=dplyr::all_of(feature_col), 
                       values_from=dplyr::all_of(value_col)) %>%
    tibble::column_to_rownames(var="specimen_id") %>%
    as.matrix()
  
  plasma_cytokine_concentration_by_legendplex_meta <-
    tibble(specimen_id=rownames(plasma_cytokine_concentration_by_legendplex_wide)) %>%
    dplyr::left_join((meta_data %>% 
                        dplyr::select(specimen_id, dataset) %>%
                        dplyr::mutate(specimen_id = as.character(specimen_id))),
                      by="specimen_id")
    
  plasma_cytokine_concentration_by_legendplex_corrected <- sva::ComBat(
    dat = t(plasma_cytokine_concentration_by_legendplex_wide),
    batch = plasma_cytokine_concentration_by_legendplex_meta$dataset
  )
  
  integrated_data$plasma_cytokine_concentration_by_legendplex <- 
    t(plasma_cytokine_concentration_by_legendplex_corrected) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="specimen_id") %>%
    tidyr::pivot_longer(cols=!c(specimen_id), 
                        names_to=feature_col,
                        values_to=value_col)
  
  # 5) plasma_cytokine_concentration_by_olink: No Normalization needed? 
  # TODO: Take into consideration that the CMI-PB consortium suggests median baseline normalization
  integrated_data$plasma_cytokine_concentration_by_olink <- 
    normalized_experimental_data$plasma_cytokine_concentration_by_olink
  
  # 6) t_cell_activation: No Normalization needed? (suggested by CMI-PB consortium)
  integrated_data$t_cell_activation <- 
    normalized_experimental_data$t_cell_activation
  
  # 7) t_cell_polarization: No Normalization needed? (suggested by CMI-PB consortium)
  integrated_data$t_cell_polarization <- 
    normalized_experimental_data$t_cell_polarization
  
  # check that we did not forget any assay
  stopifnot(dplyr::setequal(names(normalized_experimental_data), names(integrated_data)))
  
  return(integrated_data)  
}





