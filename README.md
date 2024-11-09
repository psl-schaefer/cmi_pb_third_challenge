# CMI-PB Third Challenge

- Challenge Overview: https://www.cmi-pb.org/

## Data Processing

### Clean Feature Names

- I replaced all "[+-.: /]" with other strings such that I could use the formula syntax in R

### Data Filtering

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

- Lower Limit of Quantification
  - Plasma Antibody Titers:
  	- TODO
  - Olink:
    - See forum post [here](https://discuss.cmi-pb.org/t/how-is-the-limit-of-detection-lod-estimated-for-olink-data-and-how-is-this-handled-in-the-data-analysis/122)
  	
- Outlier Removal
  - For the legendplex assay I removed these specimen: 661, 657, 658, 689, 903, 638, 659, 661
  	
- PBMC gene expression
  - For simplicity, I decided to only keep genes that have a unique mapping from ensemble id to gene symbol
  	
### Data Normalization



