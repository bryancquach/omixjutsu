# Tools for model selection prior to regression modeling.
# Author: Bryan Quach (bryancquach@gmail.com)

#' Create design matrix
#'
#' Format data table to produce a design matrix.
#'
#' Creates a design matrix and prunes variables with no variance. Removes pairwise collinear
#' variables.
#'
#' @param data A data frame of potential explanatory variables for the design matrix. Rows should
#' be records and columns should be variables.
#' @param var_names A vector of column names from `data` to consider for the design matrix.
#' @param corr_cutoff A numeric value in [0,1] denoting the absolute correlation cutoff to
#' classify pairs of collinear variables.
#' @param corr_method A string indicating the type of correlation to use for determining
#' collinearity. Must be `pearson` or `spearman`.
#' @param check_binary_vars A boolean indiciating if dichotomous variables should be included in
#' the collinearity pruning.
#' @param strings_as_factors A boolean indiciating if character variables should be converted to
#' factors in the returned data frame.
#' @return A data frame of explanatory variables for regression.
#' @seealso \code{\link{get_residual_matrix}}
#' @export
get_design_mat <- function (data,
                            var_names = NULL,
                            corr_cutoff = 0.75, 
                            corr_method = c("spearman", "pearson"),
                            check_binary_vars = T,
                            strings_as_factors = T) {
  if (corr_cutoff < 0 | corr_cutoff > 1) {
    stop("Error: 'corr_cutoff' must be in the range [0,1]")
  }
  corr_method <- match.arg(corr_method)
  if (length(var_names) > 0) {
    if (any(!var_names %in% colnames(data))) {
      stop("Error: Not all `var_names` found in `data`")
    }
    if (length(var_names) == 1) {
      data_subset <- data[, var_names, drop = F]
      if (length(unique(na.omit(data_subset[, 1, drop = T]))) > 1) {
        return(data_subset)
      } else {
        stop("Error: no `var_names` variable with non-zero variance")
      }
    } else {
      data <- data[, var_names]
    }
  }
  is_keeper <- apply(
    data, 
    2, 
    function(x) {
      num_items <- length(unique(na.omit(x)))
      return(num_items > 1)
    }
  )
  if (sum(is_keeper) > 0) {
    print(paste("Variables with non-zero variance:", sum(is_keeper)))
    if (sum(is_keeper) == 1) {
      return(data[, is_keeper, drop = F])
    }
    data <- data[, is_keeper]
    corrcheck_data <- data
  } else {
    stop("Error: no `var_names` variables with non-zero variance")
  }
  # Only binary or numeric variables can be compared for collinearity
  if (check_binary_vars) {
    is_binary <- apply(
      corrcheck_data, 
      2, 
      function(x) {
        num_items <- length(unique(na.omit(x)))
        return(num_items == 2)
      }
    )
    print(paste("Binary variables:", sum(is_binary)))
    # Binary variables must be numeric for correlation calculations
    for (i in seq_len(length(is_binary))) {
      if (is_binary[i]) {
        corrcheck_data[, i] <- as.numeric(as.factor(corrcheck_data[, i, drop = T]))
      }
    }
  } else {
    is_binary <- rep(FALSE, ncol(corrcheck_data))
  }
  is_numeric <- sapply(corrcheck_data, class) %in% c("integer", "numeric")
  corr_candidates <- which(is_numeric | is_binary)
  print(paste("Collinearity assessment candidates:", length(corr_candidates)))
  covar_mat <- Hmisc::rcorr(as.matrix(corrcheck_data[, corr_candidates]), type = corr_method)$r
  rm_variables <- caret::findCorrelation(covar_mat, cutoff = corr_cutoff, names = T)
  print("Variables to remove:")
  print(rm_variables)
  data <- data[, !colnames(data) %in% rm_variables]
  if (strings_as_factors) {
    is_character <- sapply(data, class) == "character"
    for (i in which(is_character)) {
      data[, i] <- as.factor(data[, i, drop = T])
    }
  }
  return(data)
}

#' Calculate residual matrix
#'
#' Get residuals from fitting a linear model with a matrix of outcomes.
#'
#' Calculates regression residuals using a fixed set of explanatory variables and a matrix of
#' outcomes.
#'
#' @param design A data frame with explanatory variables.
#' @param outcomes A data frame like object with outcome variables as rows and observations as
#' columns.
#' @return A list with the following items:
#' * residuals: A numeric matrix with each column representing the residuals for an outcome. The
#' ordering follows the row order of `outcomes`.
#' * dof: Model degrees of freedom.
#' @seealso \code{\link{get_matrix_r2}} \code{\link{get_design_mat}}
#' @export
get_residual_matrix <- function (design, outcomes) {
  results <- list()
  y <- outcomes
  x <- design
  y <- t(y)
  resid_mat <- matrix(ncol = ncol(y), nrow = nrow(y))
  for(i in seq_len(ncol(y))){
    lm_fit <- RcppEigen::fastLm(y[, i] ~ 1 + . , data = x)
    resid_mat[, i] <- lm_fit$residuals
  }
  results[["residuals"]] <- resid_mat
  results[["dof"]] <- lm_fit$df.residual
  return(results)
}

#' Calculate (adjusted) r-squared for a matrix
#'
#' Calculates regression model r-squared across a matrix of outcomes.
#'
#' @param outcomes A data frame like object with outcome variables as rows and observations as
#' columns.
#' @param residuals A numeric matrix with each column representing the residuals for an outcome.
#' The ordering follows the row order of `outcomes`.
#' @param adjusted A boolean indicating if the adjusted r-squared should be calculated in addition
#' to the unadjusted r-squared.
#' @param dof An integer denoting the degrees of freedom of the regression model used to compute
#' the residuals. Only required if computing adjusted r-squared.
#' @return A list with the following items:
#' * r2: A numeric value corresponding to the unadjusted r-squared.
#' * adjusted_r2: A numeric value (or NA) corresponding to the adjusted r-squared.
#' @seealso \code{\link{get_residual_matrix}}
#' @export
get_matrix_r2 <- function (outcomes,
                           residuals,
                           adjusted = F,
                           dof = NULL) {
  if (adjusted & is.null(dof)) {
    stop("Error: degrees of freedom required for adjusted r-squared calculation")
  }
  tss <- norm(outcomes - rowMeans(outcomes, na.rm = T), type = "F")^2
  rss <- norm(residuals, type = "F")^2
  n <- ncol(outcomes)
  r2 <- (tss - rss) / tss
  if (adjusted) {
    adjusted_r2 <- 1 - ((1 - r2) * ((n - 1) / dof))
  } else {
    adjusted_r2 <- NA
  }
  return(list(r2 = r2, adjusted_r2 = adjusted_r2))
}

#' Plot (adjusted) r-squared for multiple models
#'
#' Plots regression model r-squared vs. model size (number of variables).
#'
#' @param data A data frame like object with outcome variables as rows and observations as
#' columns.
#' @return A ggplot object
#' @seealso \code{\link{get_matrix_r2}} \code{\link{get_residual_matrix}} 
#' \code{\link{get_design_mat}}
#' @export
plot_r2 <- function(){}
