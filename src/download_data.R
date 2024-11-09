library(tidyverse)
library(jsonlite)
library(httr)

# Define where things should be saved
data_dir <- "data/"

# Download the demographic and experimental data
server_url <- "https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/current/3rd_challenge/"

system(paste("wget", 
             "--recursive", # download files and directories recursively
             "--no-parent", # prevents wget from ascending to parent directories
             "--no-host-directories", # prevents wget from creating a directory named after the host (e.g., www.cmi-pb.org)
             "--cut-dirs=4", # skips creating the first x levels of the directory structure from the remote URL
             "--directory-prefix", 
             "--reject 'index.html*'", # excludes any files matching index.html*, preventing wget from downloading unwanted index files
             paste0("--directory-prefix=", data_dir), # specify where to save
             server_url)
)

# Download the metadata for genes and cell types
meta_data_dir <- file.path(data_dir, "meta_data")
if (!file.exists(meta_data_dir)) {
  dir.create(meta_data_dir, recursive = TRUE)
}

external_data_dir <- file.path(data_dir, "external_data")
if (!file.exists(external_data_dir)) {
  dir.create(external_data_dir, recursive = TRUE)
}

prc_data_dir <- file.path(data_dir, "prc_datasets")
if (!file.exists(prc_data_dir)) {
  dir.create(prc_data_dir, recursive = TRUE)
}

base_url <- "https://www.cmi-pb.org/api/v5" # "https://www.cmi-pb.org/api"
endpoints <- c("/cell_type", "/gene", "/protein")
for (endpoint in endpoints) {
  url <- paste0(base_url, endpoint)
  response <- GET(url)
  if (status_code(response) == 200) {
    json_content <- content(response, as = "text", encoding = "UTF-8")
    write(json_content, file = file.path(meta_data_dir, paste0(basename(endpoint), ".json")))
    cat("Downloaded:", endpoint, "\n")
  } else {
    cat("Failed to download", endpoint, "Status code:", status_code(response), "\n")
  }
}

celltype_info <- purrr::map(jsonlite::read_json(file.path(meta_data_dir, "cell_type.json")), function(l) {
  as.data.frame(l)
}) %>%
  dplyr::bind_rows()

# NOTE: Since there is no gating info at all for 2023, 
# we just assume it is the same as in 2022
celltype_info <- dplyr::bind_rows(
  celltype_info,
  celltype_info %>%
    dplyr::filter(dataset == "2022_dataset") %>%
    dplyr::mutate(dataset = "2023_dataset")
)

# NOTE: Since there no gating info for "Basophils" and "CD3CD19neg" in the 2020 dataset, 
# we just assume it is the same as in 2021 and 2022
celltype_info <- dplyr::bind_rows(
  celltype_info,
  celltype_info %>%
    dplyr::filter(cell_type_name %in% c("Basophils", "CD3CD19neg")) %>%
    dplyr::filter(dataset == "2022_dataset") %>%
    dplyr::mutate(dataset = "2020_dataset")
)

# NOTE: Since there no gating info for "non-pDCs", 
# we just assume ..

# gating_definition_test <- c(
#   "CD19-CD3+CD14-", 
#   "CD19-CD3-CD56+/++", 
#   "CD19-CD3-CD56++"
# )
#str_split(gating_definition_test, "(?<=\\+)(?=[^+/-])|(?<=-)(?=[^+/-])")

# Clean and harmonize the gating info 
# celltype_info <- celltype_info %>%
#   dplyr::mutate(gating_definition_clean = str_replace(
#     str_replace(.data$gating_definition, "^[+\\-\\\\ ]+", ""), "[/]*([+-])[+/+-]+", "\\1")
# ) %>%
#   dplyr::mutate(
#     harmonized_gating_definition = strsplit(.data$gating_definition_clean, "(?<=\\+|-)", perl = TRUE) %>%
#       purrr::map_chr(~ paste0(sort(.x), collapse=""))
#   )
celltype_info <- celltype_info %>%
    dplyr::mutate(
      harmonized_gating_definition = str_split(.data$gating_definition, "(?<=\\+)(?=[^+/-])|(?<=-)(?=[^+/-])") %>%
        purrr::map_chr(~ paste0(sort(.x), collapse=""))
    )

# Function to compute hierarchy levels and parents based on parent-child relationships
compute_levels_and_parents <- function(harm_gating_def) {
  # gating_info <- celltype_info$gating_definition
  # gating_info <- c("CD19-CD3-", "CD19-CD3-CD56+/++", "CD19-CD3-CD56++")
  
  # Check that last character is either + or -
  stopifnot(all(purrr::map_chr(harm_gating_def, ~ str_sub(.x, -1, -1)) %in% c("+", "-")))

  # Split each string into components (split by +/-)
  split_gating_info <- strsplit(harm_gating_def, "(?<=\\+|-)", perl = TRUE)
  
  # Initialize vectors to store levels and parents
  levels <- numeric(length(harm_gating_def))
  parents <- character(length(harm_gating_def))
  
  # Loop over each string to determine its level and parent
  for (i in seq_along(harm_gating_def)) {
    # i <- 6
    # Initialize the level and parent for the current string
    current_level <- 0
    parent <- NA  # No parent by default
    
    # Compare the current string with all other gating_info
    for (j in seq_along(harm_gating_def)) {
      if (i != j) {
        # Check if string j is a parent of string i (if j's components are a prefix of i's components)
        if (all(split_gating_info[[j]] %in% split_gating_info[[i]][1:length(split_gating_info[[j]])])) {
          # Increment the level if string j is a parent
          if (levels[j] + 1 > current_level) {
            current_level <- levels[j] + 1
            parent <- harm_gating_def[j]  # Assign the parent string
          }
        }
      }
    }
    
    # Assign the calculated level and parent to the current string
    levels[i] <- current_level
    parents[i] <- parent
  }
  
  # Return a data frame with the original gating_info, their levels, and parents
  return(list(level = levels, 
              parent = parents))
}

celltype_info <- celltype_info %>%
  dplyr::group_by(dataset) %>%
  dplyr::mutate(level = compute_levels_and_parents(gating_definition)$level,
                parent_gating = compute_levels_and_parents(gating_definition)$parent) %>%
  dplyr::mutate(
    parent_name = purrr::map_chr(.data$parent_gating, ~ unique(.data$cell_type_name[.data$gating_definition==.x]))
    ) %>%
  dplyr::ungroup()

# Lastly add a column indicating whether the cell_type_name is present in the raw data
meta_df <- dplyr::bind_rows(
  read_tsv(file.path(data_dir, 
                     "harmonized_and_processed_data", 
                     "master_harmonized_data_TSV", 
                     "training_subject_specimen.tsv"), 
           show_col_types = FALSE),
  read_tsv(file.path(data_dir, 
                     "harmonized_and_processed_data", 
                     "master_harmonized_data_TSV", 
                     "challenge_subject_specimen.tsv"), 
           show_col_types = FALSE)
)
celltype_per_dataset <- list.files(path = file.path(data_dir, "raw_datasets"), 
                            recursive=TRUE, full.names=TRUE) %>%
  magrittr::extract(grepl("pbmc_cell_frequency", .)) %>%
  purrr::map(., ~ read_tsv(.x, show_col_types=FALSE)) %>%
  dplyr::bind_rows() %>%
  dplyr::left_join((meta_df %>% 
                      dplyr::select(specimen_id, dataset)),
                   by="specimen_id") %>%
  dplyr::select(dataset, cell_type_name) %>%
  dplyr::distinct()
  
celltype_info <- celltype_info %>%
  dplyr::mutate(is_in_2020 = cell_type_name %in% celltype_per_dataset$cell_type_name[celltype_per_dataset$dataset=="2020_dataset"]) %>%
  dplyr::mutate(is_in_2021 = cell_type_name %in% celltype_per_dataset$cell_type_name[celltype_per_dataset$dataset=="2021_dataset"]) %>%
  dplyr::mutate(is_in_2022 = cell_type_name %in% celltype_per_dataset$cell_type_name[celltype_per_dataset$dataset=="2022_dataset"]) %>%
  dplyr::mutate(is_in_2023 = cell_type_name %in% celltype_per_dataset$cell_type_name[celltype_per_dataset$dataset=="2022_dataset"]) %>%
  dplyr::mutate(in_experimental_data = cell_type_name %in% celltype_per_dataset$cell_type_name)

readr::write_csv2(celltype_info, file=file.path(meta_data_dir, "celltype.csv"))



gene_info <- purrr::map(jsonlite::read_json(file.path(meta_data_dir, "gene.json")), function(l) {
  tibble::tibble(
    versioned_ensembl_gene_id = l[["versioned_ensembl_gene_id"]],
    gene_symbol = l[["gene_symbol"]]
  )
}) %>%
  dplyr::bind_rows()

readr::write_csv2(gene_info, file=file.path(meta_data_dir, "gene.csv"))


protein_info <- purrr::map(jsonlite::read_json(file.path(meta_data_dir, "protein.json")), function(l) {
  tibble::tibble(
    uniprot_id = l[["uniprot_id"]],
    mapping_db = l[["mapping_db"]],
    versioned_ensembl_protein_id = l[["versioned_ensembl_protein_id"]],
    versioned_ensembl_transcript_id = l[["versioned_ensembl_transcript_id"]]
  )
}) %>%
  dplyr::bind_rows()

readr::write_csv2(protein_info, file=file.path(meta_data_dir, "protein.csv"))


# Ontology data (not sure what to do with this...)
download.file("https://www.cmi-pb.org/myterminology/cmi-pb.owl", file.path(data_dir, "cmi-pb.owl"))

library(xml2)
xml_data <- read_xml(file.path(data_dir, "cmi-pb.owl"))
xml_list <- list()

# Example: Find all 'owl:Class' elements and store their attributes and text in a list
classes <- xml_find_all(xml_data, ".//owl:Class")

# Loop through each 'Class' element and extract details
for (i in seq_along(classes)) {
  # Extract the class name or ID (could be an attribute like rdf:about or text)
  class_id <- xml_attr(classes[i], "rdf:about")  # Example using rdf:about attribute
  
  # Extract labels or other details if present (you can adjust this to your needs)
  label <- xml_text(xml_find_first(classes[i], ".//rdfs:label"))
  
  # Store the class ID and label as a named entry in the list
  xml_list[[i]] <- list(
    id = class_id,
    label = label
  )
}

# View the resulting list
#print(xml_list)


# Download the gene info often used in the recommended processing files
download.file("https://raw.githubusercontent.com/CMI-PB/second-challenge-train-dataset-preparation/main/data/gene_90_38_export.tsv", 
              file.path(meta_data_dir, "gene_90_38_export.tsv"), method = "curl")

download.file("https://zenodo.org/records/10642079/files/CMI-PB/literature_models_first_challenge-cmipb-challenge.zip?download=1",
              file.path(external_data_dir, "literature_models_first_challenge-cmipb-challenge.zip"), method="curl")

unzip(file.path(external_data_dir, "literature_models_first_challenge-cmipb-challenge.zip"), 
      exdir = file.path(external_data_dir, "literature_models_first_challenge-cmipb-challenge"))


load(file.path(external_data_dir, 
               "literature_models_first_challenge-cmipb-challenge", 
               "CMI-PB-literature_models_first_challenge-c3c23cb",
               "Study-1-Avey-2017",
               "geneSetDB.rda"))

gene_set <- strsplit(geneSetDB, "\t")      
module.names <- as.data.frame(x = sapply(gene_set,"[", 1))
colnames(module.names) <- 'module'
gene_set <- lapply(gene_set, "[",-1:-2) # remove name and description columns
gene_set <- lapply(gene_set, function(x){ x[ which( x != "") ] }) # remove empty strings
names(gene_set) <- module.names$module # adding module names to gene_set

write_rds(gene_set, file=file.path(external_data_dir, 
                                   "literature_models_first_challenge-cmipb-challenge", 
                                   "CMI-PB-literature_models_first_challenge-c3c23cb",
                                   "Study-1-Avey-2017",
                                   "gene_set.RDS"))
