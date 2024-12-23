# TODO

- Add a table of contents at the top here

# CMI-PB Third Challenge

- Challenge Overview: [https://www.cmi-pb.org/](https://www.cmi-pb.org/)

- Webpage for this project: [https://psl-schaefer.github.io/cmi_pb_third_challenge/](https://psl-schaefer.github.io/cmi_pb_third_challenge/)

## Thoughts on the Challenge

- There are several aspects that make this challenge challenging:
  - Large-P-Small-N Problem
  - Seven different prediction tasks
  - Six different data modalities
  - Strong cohort (group) effects
  
- Given the complexity of the challenge, there are many degrees of freedom of how to approach. To not get lost I tried to focus on the following aspects:
  - 1) Solid QC and data preprocessing of the experimental data.
  - 2) Solic Model evaluation framework, because it is extremly easy to overfit on the training data in the Large-P-Small-N setting.

## Data Processing

### Clean Feature Names

- I replaced all "[+-.: /]" with other strings such that I could use the formula syntax in R

### Data Filtering

- Code: `filter_experimental_data` in `read_data.R`

- Only keep Specimen for which metadata are available

- Manually subset features for certain assays:
  - PBMC frequencies: Selected a subset of cell types based on the gating info (but did not benchmark this selection): Monocytes, Classical_Monocytes, Intermediate_Monocytes, Non-Classical_Monocytes, CD8Tcells, CD4Tcells, Bcells, NK, pDC, Basophils
  - Olink: Only use proteins/cytokines that have an NA fraction below 50%.

- Olink QC Qarning: 
  - Removed all measurements with QC warning in the Olink assay (which was about 300 measurements)

- Different Units in Assay
  - Plasma Antibody Titers:
    - See forum post [here](https://discuss.cmi-pb.org/t/multiple-units-in-the-ab-titer-table/57/3)
  	- In the 2020 dataset we have `IU/ML` whereas for 2021-2023 we have MFI (fluorescence intensity). `IU/ML` is obtained using a serum standard, thus there is not trivial conversion from MFI to IU/ML.
  	- I decided not to consider the 2020 Ab titer data
  - Olink:
  	- In the 2020 dataset we have `Normalized Protein eXpression`, whereas for 2021-2023 we have `PG/ML`
  	- I removed all measurements with different units.
	- Here, I remove by far the most measurements.

- Lower Limit of Detection in Plasma Antibody Titers:
	- See forum post [here](https://discuss.cmi-pb.org/t/data-preprocessing-questions/129)
	- Removed specimen with more than 50% of measurements below LOD

- Lower Limit of Quantification in Olink assay
  - See forum post [here](https://discuss.cmi-pb.org/t/how-is-the-limit-of-detection-lod-estimated-for-olink-data-and-how-is-this-handled-in-the-data-analysis/122)
  - Removed specimen with more than 50% of measurements below LOD
  	
- Outlier Removal
  - For the legendplex assay I removed these specimen: 661, 657, 658, 689, 903, 638, 659, 661 (based on PCA plots)
  	
- PBMC gene expression
  - For simplicity, I decided to only keep genes that have a unique mapping from ensemble id to gene symbol
  	
### Data Normalization

- Code: `normalize_experimental_data` in `normalize_integrate.R`

- pbmc_cell_frequency: Median Baseline Normalization (i.e. divide each feature by the median of that feature in the measurments from specimen from day 0)

- pbmc_gene_expression: VST

- plasma_ab_titer: Median Baseline Normalization

- plasma_cytokine_concentration_by_legendplex: Median Baseline Normalization

- plasma_cytokine_concentration_by_olink: No Normalization

- t_cell_activation: No Normalization

- t_cell_polarization: No Normalization

### Data Integration

- Code: `integrate_experimental_data` in `normalize_integrate.R`

- pbmc_cell_frequency: No Integration

- pbmc_gene_expression: ComBat-seq

- plasma_ab_titer: No Integration

- plasma_cytokine_concentration_by_legendplex: ComBat

- plasma_cytokine_concentration_by_olink: No Integration

- t_cell_activation: No Integration

- t_cell_polarization: No Integration

### Further Processing

- pbmc_gene_expression: Selection of the `2000` most highly variable genes (simply using dispersion)

- Missing values were imputed per assay using `impute::impute.knn` with `k=7`

## Construction of Target Tracks

- Code: `generate_all_targets` in `generate_targets.R`

- Which subjects / specimens did I consider:

- I defined the following acceptable differences between the actual and planned day of visit wrt booster administration:
  - day_0 = list("min_diff" = -7, "max_diff" = 0)
  - day_1 = list("min_diff" = 1, "max_diff" = 2)
  - day_3 = list("min_diff" = 3, "max_diff" = 6)
  - day_14 = list("min_diff" = 12, "max_diff" = 16)
  - day_30 = list("min_diff" = 28, "max_diff" = 40)

- In tasks 1.2, 2.2, 3.2, instead of predicting fold changes, I predicted log fold changes, as the distributions of the logFC looked less skewed. 

- In task 3.1, instead of predicting the raw expression of CCL3, I predicted the log expression

## Model Selection

- Code: `model_selection.qmd`

- Given that we deal here with a Large-P-Small-N problem, having a good model selection framework is probably even more important than the actual model.

- To evaluate and select models for each task, I used nested cross validation (CV):
  - Folds in the outer loop: Each cohort is a fold (i.e. Group k-fold) -> Outer loop is used to estimate cross-cohort model performance
  - Fold in the inner loop: Each subject is a fold (i.e. LOOCV) -> Inner loop is used to select the best set of hyperparameters for any model. Using this set on hyperparameters I then estimate model performance on the hold-out fold from the outer loop

- I used this model evaluation framework to test many combinations of models and features:
  - Models: LASSO, Elastic Net, Random Forest
    - XGboost (or any other kind of boosting algorithm) is usually state-of-the-art in tabular ML problems, but here I was not sure about using it given the small amount of training data points we have (large-p-small-n), ...
  - Features: Power set of all assays

- To choose the final model per task, I not only considered the perfomrance, but all of the following points:
  - Top ranking performance
  - Low variance between the test sets
  - If several models with top performance and low variance, choose the more regularized model
  - If several models with top performance and low variance, choose model that is not using any assay that is missing often in test data (specifically Olink!)

## How to Run the Code

- Clone/Fork the repo, e.g. `git clone https://github.com/psl-schaefer/cmi_pb_third_challenge`

- Cd into the cloned repo (using your favorite shell)

- Run `Rscript -e "renv::restore()"`

- Check with `Rscript -e "renv::status()"`
  - Should return: `No issues found -- the project is in a consistent state.`
  
- Alternatively, if `renv` does not work, one could also delete the lock file and renv directory, and manually install the packages.

- Run `Rscript src/download_data.R`

- Run `Rscript src/compute_gene_stats.R`

- a) Generating the whole webpage (i.e. all analyses):
  - Install quarto
  - Note, if you want to include the MOFA model (which I ended up not using), change the specification for the conda environment that is used by `reticulate`. E.g. whereever you see this line `reticulate::use_condaenv("scanpy")`, use a conda environment that is installed in your system.
  - Run `quarto render` inside the project directory (in your favorite shell)
  
- b) Generating only the predictions
  - Run `quarto render model_selection.qmd` (on my machine this runs about 4 hours)
  
- Submission file should be in `results/model_selection_default/` (saved using the current data as prefix)




