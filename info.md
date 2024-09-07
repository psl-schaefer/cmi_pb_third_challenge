
# Overview

- Description: https://www.cmi-pb.org/blog/prediction-challenge-overview/

# Data and resources

## Study design and multi-omics datasets

Our cohort is composed of acellular-Pertussis (aP) vs whole-cell Pertussis (wP) infancy-primed subjects boosted with Tdap. We recruited individuals born before 1995 (wP) and after 1996 (aP), collected baseline blood samples, provided the Tdap booster vaccine, and then obtained blood at 1, 3, 7, 14, and 28 days post booster vaccination. We recommend reviewing the 'Learn about the project' and 'Understand the data' sections to delve deeper into study design, experimental data generation and standardization.

From the samples, we generated omics data by:

Cell frequencies of PBMCs (30+ cell populations) using flow cytometry,
PBMC gene expression (50,000+ genes),
Plasma cytokine and chemokine concentrations (45+ proteins) using Olink and LegendPlex
Plasma antigen-specific antibody measurements (multiple antigens, total IgG, and subclasses of IgG)
T cell polarization using FluoroSpot assay
T cell activation using AIM assay

Challenge data

The data has been split into two groups:

Training dataset (Baseline and longitudinal readouts: 2020, 2021, and 2022 dataset). To build your computational models, use the training set, which includes the outcome (also known as the "ground truth") for each subject. Your model will be based on features extracted from longitudinal multi-omics readouts (seven assays) as well as demographic data (such as age, infancy vaccination, and biological sex). You can create new features using feature engineering techniques.

Challenge dataset (Baseline readouts: 2023 dataset). Use the challenge dataset to evaluate your model's performance on new, unseen baseline data. The challenge dataset does not provide the ground truth (longitudinal vaccine response) for each subject. Your task is to predict these outcomes using the model you built. Use your model to predict the vaccine response for each subject in the challenge dataset.

- Link: https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/current/

```bash
wget -r -np -nH --cut-dirs=2 -R "index.html*" https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/current/
mkdir data
mv -r current/3rd_challenge/* data/
```

