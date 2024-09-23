library(tidyverse)
library(ggdark)
library(vsn)
library(lme4)
library(matrixStats)

pvca_analysis <- function(assay_mtx_wide, meta_data) {
  
  # center the data, prior to pca
  assay_mtx_wide_centered <- 
    t((t(assay_mtx_wide) - colMeans(assay_mtx_wide)))
  #stopifnot(all(dplyr::near(colMeans(assay_mtx_wide_centered), 0)))
  if (!all(dplyr::near(colMeans(assay_mtx_wide_centered), 0))) {
    warning("Features are not perfectly centered due to ...")
  }
  # assay_mtx_wide_centered_scaled <-
  #   t(t(assay_mtx_wide_centered) / matrixStats::colSds(assay_mtx_wide_centered))
  # stopifnot(all(dplyr::near(matrixStats::colSds(assay_mtx_wide_centered_scaled), 1)))
  
  # compute full PCA
  #S <- 1/nrow(assay_mtx_wide_centered_scaled) * t(assay_mtx_wide_centered_scaled) %*% assay_mtx_wide_centered_scaled
  theDataCor <- cor(t(assay_mtx_wide_centered)) # why do we compute the covariance between instances?
  eigenData <- eigen(theDataCor)
  eigenValues = eigenData$values
  ev_n <- length(eigenValues)
  eigenVectorsMatrix = eigenData$vectors
  eigenValuesSum = sum(eigenValues)
  percents_PCs = eigenValues/eigenValuesSum
  
  # get the covariates to regress against pca components
  exp_design <- meta_data %>%
    dplyr::select(specimen_id, infancy_vac, biological_sex, ethnicity, age_at_boost, dataset) %>%
    dplyr::distinct() %>%
    dplyr::mutate(specimen_id = as.character(specimen_id)) %>%
    tibble::column_to_rownames(var="specimen_id") %>%
    as.data.frame() %>%
    .[rownames(assay_mtx_wide_centered), ]
  
  expDesignRowN <- nrow(exp_design)
  expDesignColN <- ncol(exp_design)
  
  threshold <- 0.8
  n_pcs = sum(cumsum(percents_PCs) < threshold)
  n_pcs <- ifelse(n_pcs < 3, 3, n_pcs)
  
  # flatten the first n_pcs into the pc_data_matrix (which is just a column vector)
  pc_data_matrix <- matrix(eigenVectorsMatrix[, 1:n_pcs], 
                           nrow = expDesignRowN * n_pcs, 
                           ncol = 1)
  
  # why are we repeating this?
  AAA <- exp_design[rep(1:expDesignRowN, n_pcs), ]
  
  Data <- cbind(AAA, pc_data_matrix)
  
  variables <- colnames(exp_design)
  
  for (i in 1:length(variables)) {
    Data$variables[i] <- as.factor(Data$variables[i])
  }
  
  op <- options(warn = (-1))
  effects_n = expDesignColN + choose(expDesignColN, 2) + 1
  
  randomEffectsMatrix <- matrix(data = 0, nrow = n_pcs, ncol = effects_n)
  model.func <- c()
  index <- 1
  for (i in 1:length(variables)) {
    mod = paste("(1|", variables[i], ")", sep = "")
    model.func[index] = mod
    index = index + 1
  }
  for (i in 1:(length(variables) - 1)) {
    for (j in (i + 1):length(variables)) {
      mod = paste("(1|", variables[i], ":", variables[j], 
                  ")", sep = "")
      model.func[index] = mod
      index = index + 1
    }
  }
  function.mods <- paste(model.func, collapse = " + ")
  for (i in 1:n_pcs) {
    y = (((i - 1) * expDesignRowN) + 1)
    funct <- paste("pc_data_matrix", function.mods, sep = " ~ ")
    Rm1ML <- lmer(funct, Data[y:(((i - 1) * expDesignRowN) + 
                                   expDesignRowN), ], REML = TRUE, verbose = FALSE, 
                  na.action = na.omit)
    randomEffects <- Rm1ML
    randomEffectsMatrix[i, ] <- c(unlist(VarCorr(Rm1ML)), 
                                  resid = sigma(Rm1ML)^2)
  }
  effectsNames <- c(names(getME(Rm1ML, "cnms")), "resid")
  randomEffectsMatrixStdze <- matrix(data = 0, nrow = n_pcs, 
                                     ncol = effects_n)
  for (i in 1:n_pcs) {
    mySum = sum(randomEffectsMatrix[i, ])
    for (j in 1:effects_n) {
      randomEffectsMatrixStdze[i, j] = randomEffectsMatrix[i, 
                                                           j]/mySum
    }
  }
  randomEffectsMatrixWtProp <- matrix(data = 0, nrow = n_pcs, 
                                      ncol = effects_n)
  for (i in 1:n_pcs) {
    weight = eigenValues[i]/eigenValuesSum
    for (j in 1:effects_n) {
      randomEffectsMatrixWtProp[i, j] = randomEffectsMatrixStdze[i, 
                                                                 j] * weight
    }
  }
  randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
  randomEffectsSums <- colSums(randomEffectsMatrixWtProp)
  totalSum = sum(randomEffectsSums)
  randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, 
                                         ncol = effects_n)
  for (j in 1:effects_n) {
    randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum
  }
  var_expl_df <- tibble(variance_expl = as.vector(randomEffectsMatrixWtAveProp), 
                        covariate = factor(effectsNames)) %>%
    dplyr::arrange(desc(variance_expl))
  
  p1 <- var_expl_df %>%
    ggplot(mapping = aes(y = reorder(covariate, variance_expl), x = variance_expl)) +
    geom_bar(stat = "identity", fill="blue") +
    geom_text(aes(label = signif(variance_expl, digits = 2)),
              nudge_y   = 0.01,
              size      = 4) +
    labs(y = "Covariate", y = "Proportion of the variance explained") + 
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = 1,
                                     size  = 12)) +
    ggdark::dark_mode(verbose=FALSE)
  list(var_expl_df = var_expl_df, 
       p1 = p1)
}


generate_pca_plots <- function(assay_mtx_wide, meta_data) {
  # assay_mtx_wide <- wide_experimental_data$pbmc_cell_frequency
  set.seed(42)
  pca_out <- prcomp(assay_mtx_wide, center=TRUE, retx=TRUE, scale.=TRUE)
  
  rownames(pca_out$rotation) <- colnames(assay_mtx_wide)
  rownames(pca_out$x) <- rownames(assay_mtx_wide)
  
  assay_mtx_wide_df <- as.data.frame(assay_mtx_wide) %>%
    tibble::rownames_to_column(var="specimen_id") %>%
    dplyr::mutate(specimen_id=as.integer(specimen_id))
  
  dataset_fill <- assay_mtx_wide_df %>%
    dplyr::left_join(meta_data, by="specimen_id") %>%
    dplyr::pull(dataset)
  
  plot_df <- pca_out$x %>%
    as.data.frame() %>%
    rownames_to_column("specimen_id") %>%
    dplyr::mutate(specimen_id=as.integer(specimen_id)) %>%
    dplyr::left_join(assay_mtx_wide_df, by="specimen_id") %>%
    dplyr::left_join(meta_data, by="specimen_id")
  
  var_contrib <- get_pca_var(pca_out)$contrib
  top_variables <- rowSums(var_contrib) %>%
    sort(decreasing=TRUE) %>%
    head(10) %>%
    names()
  
  p1 <- plot_df %>%
    ggplot() +
    geom_point(aes(x=PC1, y=PC2, color=dataset),
               alpha=0.5) + 
    coord_equal() +
    ggdark::dark_mode(verbose=FALSE)
  
  p2 <- plot_df %>%
    dplyr::select(dplyr::all_of(c("PC1", "PC2", top_variables))) %>%
    tidyr::pivot_longer(cols=dplyr::all_of(top_variables), 
                        values_to = "values",
                        names_to = "features") %>%
    ggplot() +
    geom_point(aes(x=PC1, y=PC2, color=values), alpha=0.5) + 
    facet_wrap(~features) +
    coord_equal() +
    ggdark::dark_mode(verbose=FALSE) +
    scale_color_viridis_c()
  
  p3 <- fviz_pca_biplot(pca_out, 
                        geom.ind = "point",
                        pointshape = 21,
                        pointsize = 2.5,
                        fill.ind = dataset_fill,
                        col.ind = "black",
                        repel = TRUE,
                        select.var = list(name = top_variables)
  ) +
    coord_equal() +
    ggdark::dark_mode(verbose=FALSE)
  return(
    list(p1=p1, p2=p2, p3=p3)
  )
}
