
experimental_data_settings <- list(
  "pbmc_cell_frequency" = list(
    "feature_col" = "cell_type_name",
    "value_col" = "percent_live_cell",
    "feature_subset" = c("Monocytes", # no parent, three children: Classical_Monocytes, Intermediate_Monocytes, Non-Classical_Monocytes
                         "Classical_Monocytes", # only parent: Monocytes; no children
                         "Intermediate_Monocytes", # only parent: Monocytes; no children
                         "Non-Classical_Monocytes", # only parent: Monocytes; no children
                         "CD8Tcells", # four children: NaiveCD8, TcmCD8, TemCD8, TemraCD8
                         "CD4Tcells", # four children: NaiveCD4, TcmCD4, TemCD4, TemraCD4
                         "Bcells", # only has one child: B cells (CD19+CD3-CD14-CD56-)
                         "NK", # only has one child: CD56high NK cells
                         "pDC", # has no children and no parent
                         "Basophils", # has no children and no parent
                         "CD56+CD3+T cells" # no children; parent: T cells (CD19-CD3+CD14-) (which is not present in the data)
                         )
  ),
  "pbmc_gene_expression" = list(
    "feature_col" = "versioned_ensembl_gene_id",
    "value_col" = "raw_counts"
  ),
  "plasma_antibody_levels" = list(
    "feature_col" = "isotype_antigen",
    "value_col" = "MFI_normalised"
  ),
  "plasma_cytokine_concentrations_by_legendplex" = list(
    "feature_col" = "protein_id",
    "value_col" = "concentration"
  ),
  "plasma_cytokine_concentrations_by_olink" = list(
    "feature_col" = "protein_id",
    "value_col" = "concentration"
  ),
  "t_cell_activation" = list(
    "feature_col" = "stimulation",
    "value_col" = "analyte_percentages"
  ),
  "t_cell_polarization" = list(
    "feature_col" = "protein_id",
    "value_col" = "analyte_counts"
  )
)

# celltype_info %>% 
#   dplyr::filter(level %in% c(0,1, 2)) %>% 
#   dplyr::filter(in_experimental_data) %>%
#   dplyr::select(cell_type_name, level, gating_definition, harmonized_gating_definition, parent_gating, parent_name, in_experimental_data) %>%
#   dplyr::distinct() %>%
#   View()

# pbmc_data <- dplyr::bind_rows(
#   read_tsv(file.path(data_dir, 
#                      "harmonized_and_processed_data", 
#                      "master_harmonized_data_TSV", 
#                      "training_pbmc_cell_frequency_long.tsv"), 
#            show_col_types = FALSE),
#   read_tsv(file.path(data_dir, 
#                      "harmonized_and_processed_data", 
#                      "master_harmonized_data_TSV", 
#                      "challenge_pbmc_cell_frequency_long.tsv"), 
#            show_col_types = FALSE)
# ) %>%
#   dplyr::left_join((meta_df %>% 
#                       dplyr::select(specimen_id, dataset)),
#                    by="specimen_id")
