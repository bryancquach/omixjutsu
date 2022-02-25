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
#' Apply M-value transform to DNA methylation (DNAm) array beta intensity values.
#'
#' @param data A numeric vector or matrix of DNAm array beta intensity values.
#' @return A numeric vector or matrix of M-values.
#' @export
beta_to_mvalue <- function(data) {
  mvalues <- log2(data / (1 - data))
  return(mvalues)
}

#' DNA methylation beta values
#'
#' Convert DNA methylation (DNAm) array M-values to beta intensity values.
#'
#' @param data A numeric vector or matrix of DNAm array M-values.
#' @return A numeric vector or matrix of beta intensity values.
#' @export
mvalue_to_beta <- function(data) {
  betas <- (2 ^ data) / ((2 ^ data) + 1)
  return(betas)
}
