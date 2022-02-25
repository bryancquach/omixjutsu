# Functions for applying transformations to data vectors or matrices.
# Author: Bryan Quach (bryancquach@gmail.com)

#' Rank inverse-normal transform
#'
#' Apply the rank inverse-normal transform (RINT) to a vector.
#'
#' @param data A numeric vector.
#' @return A numeric vector of normalized values.
#' @export
rint <- function(data) {
  if (!is.vector(data)) {
    stop("Error: Input data is not a vector.")
  }
  zscore <- scale(data)
  rankings <- rank(zscore) - 0.5
  scaled_rank <- rankings / (max(rankings, na.rm = T) + 0.5)
  rint_values <- qnorm(scaled_rank)
  return(rint_values)
}

#' DNA methylation M-values
#'
#' Apply M-value transform to DNA methylation (DNAm) array beta intensity values. Provides range
#' limit options to prevent calculation of -Inf and Inf values.
#'
#' @param data A numeric vector or matrix of DNAm array beta intensity values.
#' @param min_value An optional numeric value denoting the min M-value allowed. Any M-values that
#' are less than this value will be set to this value.
#' @param max_value An optional numeric value denoting the max M-value allowed. Any M-values that
#' exceed this value will be set to this value.
#' @return A numeric vector or matrix of M-values.
#' @export
beta_to_mvalue <- function(data, min_value = -Inf, max_value = Inf) {
  mvalues <- log2(data / (1 - data))
  min_outliers <- which(mvalues < min_value)
  max_outliers <- which(mvalues > max_value)
  if (length(min_outliers) > 0) {
    mvalues[min_outliers] <- min_value
  }
  if (length(max_outliers) > 0) {
    mvalues[max_outliers] <- max_value
  }
  return(mvalues)
}

#' DNA methylation beta values
#'
#' Convert DNA methylation (DNAm) array M-values to beta intensity values.
#'
#' @param data A numeric vector or matrix of DNAm array M-values.
#' @param min_value An optional numeric value denoting the min M-value threshold. Any M-values that
#' are equal or less than this will have an assigned beta value of 0.
#' @param max_value An optional numeric value denoting the max M-value threhsold. Any M-values that
#' are equal or exceed this value will have an assigned beta value of 1.
#' @return A numeric vector or matrix of beta intensity values.
#' @export
mvalue_to_beta <- function(data, min_value = -Inf, max_value = Inf) {
  inverse_log_values <- 2 ^ data
  min_outliers <- which(inverse_log_values <= 2 ^ min_value)
  max_outliers <- which(inverse_log_values >= 2 ^ max_value)
  betas <- (inverse_log_values) / ((inverse_log_values) + 1)
  if (length(min_outliers) > 0) {
    betas[min_outliers] <- 0
  }
  if (length(max_outliers) > 0) {
    betas[max_outliers] <- 1
  }
  return(betas)
}
