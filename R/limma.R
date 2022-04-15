# Tools for differential gene expression (DGE) analysis using limma.
# Author: Bryan Quach (bryancquach@gmail.com)

#' @rdname low_expression_filter
#' @return A filtered DGEList object if a DGEList object is provided as input.
#' @export
low_expression_filter.DGEList <- function(dge_data,
                                          value_cutoff,
                                          min_sample_fraction = NULL,
                                          threshold_variable = NULL,
                                          use_min_fraction = T,
                                          approximate = F) {
  if (is.null(min_sample_fraction) & is.null(threshold_variable)) {
    stop("Error: either 'min_sample_fraction' or 'threshold_variable' must not be NULL")
  }
  count_matrix <- dge_data$counts
  if (is.null(min_sample_fraction)) {
    pheno_data <- dge_data$samples
    variable_freq <- table(pheno_data[, threshold_variable])
    variable_fraction <- variable_freq / nrow(pheno_data)
    if (use_min_fraction) {
      sample_cutoff <- min(variable_fraction)
    } else {
      sample_cutoff <- max(variable_fraction)
    }
  } else {
    sample_cutoff <- min_sample_fraction
  }
  filtered_count_matrix <- low_expression_filter.matrix(
    object = count_matrix,
    value_cutoff = value_cutoff,
    sample_cutoff = sample_cutoff,
    approximate = approximate
  )
  return(dge_data[rownames(filtered_count_matrix), ])
}

#' @title Get counts from DGEList object
#' @param dge_data A DGEList object.
#' @param normalized A boolean denoting if counts should be normalized.
#' @return A numeric matrix.
#' @export
counts.DGEList <- function(object, normalized = F) {
  if (normalized) {
    if (all(object$samples$norm.factors == 1)) {
      warning("All norm factors are 1. Have you called `calcNormFactors` on this DGEList object?")
    }
    scaling_factor <- object$samples$lib.size * object$samples$norm.factors
    norm_counts <- t(t(object$counts) / scaling_factor) * mean(scaling_factor)
    return(norm_counts)
  } else {
    return(object$counts)
  }
}

#' @importFrom BiocGenerics counts
#' @export
setMethod("counts", signature(object = "DGEList"), counts.DGEList)

#' @rdname per_sample_count_distribution
#' @export
per_sample_count_distribution.DGEList <- function(dge_data,
                                                  normalized = T,
                                                  point_size = 2.5,
                                                  point_alpha = 1,
                                                  y_lim = NULL) {
  if (normalized) {
    scaling_factor <- dge_data$samples$lib.size * dge_data$samples$norm.factors
    norm_counts <- t(t(dge_data$counts) / scaling_factor) * mean(scaling_factor)
    log_counts <- log10(norm_counts + 1)
  } else {
    log_counts <- log10(dge_data$counts + 1)
  }
  output_plot <- per_sample_count_distribution.matrix(
    data = log_counts,
    point_size = point_size,
    point_alpha = point_alpha,
    y_lim = y_lim
  )
  return(output_plot)
}
