library(tidyverse)
library(DESeq2)

normalize_experimental_data <- function(meta_data, 
                                        experimental_data, 
                                        gene_meta) {
  specimen_per_day <- get_specimen_per_day(meta_data=meta_data)
  normalized_data <- list()
  
  # 1) pbmc_cell_frequency: Median Baseline Normalization
  baseline_medians <- experimental_data$pbmc_cell_frequency %>%
    dplyr::inner_join(specimen_per_day$day_0, by="specimen_id") %>%
    dplyr::group_by(dataset, cell_type_name) %>%
    dplyr::summarise(median = median(percent_live_cell, na.rm=TRUE), 
                     .groups="drop") %>%
    dplyr::mutate(median = ifelse(dplyr::near(median, 0), 1, median)) # prevent division by 0?
  
  normalized_data$pbmc_cell_frequency <- experimental_data$pbmc_cell_frequency %>%
    dplyr::left_join((meta_data %>% dplyr::select(specimen_id, dataset)),
                     by="specimen_id") %>%
    dplyr::left_join(baseline_medians, by=c("dataset", "cell_type_name")) %>%
    dplyr::mutate(percent_live_cell = percent_live_cell / median)
  
  # 2) pbmc_gene_expression: VST
  # a) filter for meaningful genes
  feature_col <- experimental_data_settings$pbmc_gene_expression$feature_col
  value_col <- experimental_data_settings$pbmc_gene_expression$value_col
  pbmc_gene_expression_wide <- experimental_data$pbmc_gene_expression %>%
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
  
  # only keep genes that have a unique mapping from ensemble id to gene symbol
  ensemble_to_gene <- tibble(versioned_ensembl_gene_id = colnames(pbmc_gene_expression_wide)) %>%
    dplyr::left_join(gene_meta, by="versioned_ensembl_gene_id") %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::arrange(desc(n), gene_symbol) %>%
    dplyr::filter(n==1) %>%
    dplyr::ungroup()
  pbmc_gene_expression_wide <- pbmc_gene_expression_wide[ , ensemble_to_gene$versioned_ensembl_gene_id]
  colnames(pbmc_gene_expression_wide) <- ensemble_to_gene$gene_symbol
  
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
  pbmc_gene_expression_normalized_long <- pbmc_gene_expression_normalized %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="specimen_id") %>%
    tidyr::pivot_longer(cols=!c(specimen_id), 
                        names_to=feature_col,
                        values_to=value_col)
  normalized_data$pbmc_gene_expression <- 
    pbmc_gene_expression_normalized_long
  
  # 3) plasma_ab_titer: Median Baseline Normalization
  baseline_medians <- experimental_data$plasma_ab_titer %>%
    dplyr::inner_join(specimen_per_day$day_0, by="specimen_id") %>%
    dplyr::group_by(dataset, isotype_antigen) %>%
    dplyr::summarise(median = median(MFI, na.rm=TRUE), .groups="drop") %>%
    dplyr::mutate(median = ifelse(dplyr::near(median, 0), 1, median)) # prevent division by 0?
  
  normalized_data$plasma_ab_titer <- experimental_data$plasma_ab_titer %>%
    dplyr::left_join((meta_data %>% dplyr::select(specimen_id, dataset)),
                     by="specimen_id") %>%
    dplyr::left_join(baseline_medians, by=c("dataset", "isotype_antigen")) %>%
    dplyr::mutate(MFI = MFI / median)
  
  # 4) plasma_cytokine_concentration_by_legendplex: Median Baseline Normalization
  baseline_medians <- experimental_data$plasma_cytokine_concentration_by_legendplex %>%
    dplyr::inner_join(specimen_per_day$day_0, by="specimen_id") %>%
    dplyr::group_by(dataset, protein_id) %>%
    dplyr::summarise(median = median(concentration, na.rm=TRUE), 
                     .groups="drop") %>%
    dplyr::mutate(median = ifelse(dplyr::near(median, 0), 1, median)) # prevent division by 0?
  
  normalized_data$plasma_cytokine_concentration_by_legendplex <- 
    experimental_data$plasma_cytokine_concentration_by_legendplex %>%
    dplyr::left_join((meta_data %>% dplyr::select(specimen_id, dataset)),
                     by="specimen_id") %>%
    dplyr::left_join(baseline_medians, by=c("dataset", "protein_id")) %>%
    dplyr::mutate(concentration = concentration / median)
  
  # 5) plasma_cytokine_concentration_by_olink: No Normalization needed? 
  # TODO: Take into consideration that the CMI-PB consortium suggests median baseline normalization
  normalized_data$plasma_cytokine_concentration_by_olink <- 
    experimental_data$plasma_cytokine_concentration_by_olink
  
  # 6) t_cell_activation: No Normalization needed? (suggested by CMI-PB consortium)
  normalized_data$t_cell_activation <- 
    experimental_data$t_cell_activation
  
  # 7) t_cell_polarization: No Normalization needed? (suggested by CMI-PB consortium)
  normalized_data$t_cell_polarization <- 
    experimental_data$t_cell_polarization
  
  # check that we did not forget any assay
  stopifnot(dplyr::setequal(names(experimental_data), names(normalized_data)))
  
  return(normalized_data)
}
