# Functions for applying transformations to data vectors.
# Author: Bryan Quach (bryancquach@gmail.com)

#' Rank inverse-normal transform
#'
#' Apply the rank inverse-normal transform (RINT) to a vector.
#'
#' @param data A numeric vector.
#' @return A numeric vector.
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
