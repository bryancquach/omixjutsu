# Convenience functions for working with tabular data
# Author: Bryan Quach (bryancquach@gmail.com)

#' Exclude table dimenions by matching values.
#'
#' Subset table to exclude rows/columns that match a specific value.
#'
#' Uses a user-specified row or column to find matches to a user-specified value. For matches on a
#' row, corresponding columns are removed. For matches on a column, corresponding rows are removed.
#'
#' @param data Data frame to subset.
#' @param dimension The dimension from which to pull values for matching. '1' for rows and '2' for
#' columns. The dimension not specified is the one that will be reduced in length.
#' @param by The row/column name or index to use to pull values for matching.
#' @param pattern The pattern to match against. Can be any fixed string, `NA`, or `NULL` values.
#' @param trim_whitespace If `TRUE`, remove leading and trailing whitespace from values to match.
#' @return A data frame with one of the dimensions reduced to exclude rows/columns that match
#' against the user-specified value.
#' @export
filter_table <- function(data, dimension, by, pattern, trim_whitespace = T) {
  if (!dimension %in% c(1, 2)) {
    stop("Error: 'dimension' must be '1' (row) or '2' (column)")
  }
  index <- by
  if (is.character(by)) {
    if (dimension == 1) {
      index <- which(rownames(data) == by)
    } else {
      index <- which(colnames(data) == by)
    }
  }
  if (dimension == 1) {
    data_vec <- unlist(data[index, ])
  } else {
    data_vec <- unlist(data[, index])
  }
  if (is.factor(data_vec)) {
    data_vec <- as.character(data_vec)
  }
  if (trim_whitespace) {
    data_vec <- trimws(data_vec)
  }
  if (is.null(pattern)) {
    rm_index <- which(is.null(data_vec))
  } else if (is.na(pattern)) {
    rm_index <- which(is.na(data_vec))
  } else {
    rm_index <- which(data_vec == pattern)
  }
  # Early return when no matches are found
  if (length(rm_index) < 1) {
    return(data)
  }
  if (dimension == 1) {
    data_subset <- data[, -1 * rm_index, drop = F]
  } else {
    data_subset <- data[-1 * rm_index, , drop = F]
  }
  return(data_subset)
}

#' Detect samples above a specified cutoff value
#'
#' Identify samples with values that exceed a minimum cutoff value.
#'
#' Indicates whether a sample exceeds a minimum cutoff value for a user-specified metric. Includes
#' options to set the cutoff using values derived from interquartile range for both skewed and
#' unskewed distributions.
#'
#' @param metrics_table Data frame or matrix with metrics as columns and samples as rows.
#' @param metric Name of column with metric to use.
#' @param cutoff Cutoff value.
#' @param include_equals If `TRUE`, designate metric values equal to cutoff value as TRUE.
#' @param na_as_false If `TRUE`, designate comparisons to 'NA' values as FALSE.
#' @param coef Constant to multiply the interquartile range by in calculating a relative cutoff.
#' The default value of 1.7 approximately equates to values at least 3 standard deviations from
#' the mean for a gaussian distribution.
#' @param use_unskewed_iqr If `TRUE`, set a relative cutoff value using the interquartile range
#' of the metric distribution multiplied by `coef`. More suitable for non-skewed distributions.
#' @param use_skewed_iqr If `TRUE`, set relative cutoff value using a medcouple-based adjustment
#' of the interquartile range of the metric distribution multipled by `coef`. More suitable for
#' skewed distributions.
#' @return A list containing the following elements:
#' * `keep` a logical vector indicating whether a given value for a metric meets the designated
#' cutoff criteria. The vector order matches the input row order of `metrics_table`.
#' * `cutoff` The cutoff value used.
#' @export
exceeds_min_cutoff <- function(metrics_table,
                               metric,
                               cutoff = NULL,
                               include_equals = F,
                               na_as_false = F,
                               coef = 1.7,
                               use_unskewed_iqr = F,
                               use_skewed_iqr = F) {
  if (is.null(cutoff) & (!use_unskewed_iqr & !use_skewed_iqr)) {
    stop("Error: neither absolute nor relative cutoff value specified")
  }
  if (!is.null(cutoff) & (use_unskewed_iqr | use_skewed_iqr)) {
    stop("Error: absolute and relative cutoff both specified")
  }
  if (use_unskewed_iqr & use_skewed_iqr) {
    stop("Error: both types of relative cutoffs specified")
  }
  if (!metric %in% colnames(metrics_table)) {
    stop("Error: metric name not found in 'metrics_table'")
  }
  output <- list(keep = NA, cutoff = cutoff)
  metric_vec <- as.numeric(unlist(metrics_table[, metric]))
  if (use_unskewed_iqr) {
    output$cutoff <- quantile(metric_vec, probs = 0.25, na.rm = T) - coef * IQR(metric_vec)
  }
  if (use_skewed_iqr) {
    adjbox_stats <- robustbase::adjboxStats(
      metric_vec,
      coef = coef,
      do.conf = F,
      do.out = F
    )
    output$cutoff <- adjbox_stats$fence[1]
  }
  if (include_equals) {
    output$keep <- metric_vec >= output$cutoff
  } else {
    output$keep <- metric_vec > output$cutoff
  }
  if (na_as_false) {
    if (sum(is.na(output$keep)) > 0) {
      output$keep[is.na(output$keep)] <- F
    }
  }
  return(output)
}

#' Detect samples below a specified cutoff value
#'
#' Identify samples with values below a maximum cutoff value.
#'
#' Indicates whether a sample is below a maximum cutoff value for a user-specified metric. Includes
#' options to set the cutoff using values derived from interquartile range for both skewed and
#' unskewed distributions.
#'
#' @inheritParams exceeds_min_cutoff
#' @return A list containing the following elements:
#' * `keep` a logical vector indicating whether a given value for a metric meets the designated
#' cutoff criteria. The vector order matches the input row order of `metrics_table`.
#' * `cutoff` The cutoff value used.
#' @export
below_max_cutoff <- function(metrics_table,
                             metric,
                             cutoff = NULL,
                             include_equals = F,
                             na_as_false = F,
                             coef = 1.7,
                             use_unskewed_iqr = F,
                             use_skewed_iqr = F) {
  if (is.null(cutoff) & (!use_unskewed_iqr & !use_skewed_iqr)) {
    stop("Error: neither absolute nor relative cutoff value specified")
  }
  if (!is.null(cutoff) & (use_unskewed_iqr | use_skewed_iqr)) {
    stop("Error: absolute and relative cutoff both specified")
  }
  if (use_unskewed_iqr & use_skewed_iqr) {
    stop("Error: both types of relative cutoffs specified")
  }
  if (!metric %in% colnames(metrics_table)) {
    stop("Error: metric name not found in 'metrics_table'")
  }
  output <- list(keep = NA, cutoff = cutoff)
  metric_vec <- as.numeric(unlist(metrics_table[, metric]))
  if (use_unskewed_iqr) {
    output$cutoff <- quantile(metric_vec, probs = 0.75, na.rm = T) + coef * IQR(metric_vec)
  }
  if (use_skewed_iqr) {
    adjbox_stats <- robustbase::adjboxStats(
      metric_vec,
      coef = coef,
      do.conf = F,
      do.out = F
    )
    output$cutoff <- adjbox_stats$fence[2]
  }
  if (include_equals) {
    output$keep <- metric_vec <= output$cutoff
  } else {
    output$keep <- metric_vec < output$cutoff
  }
  if (na_as_false) {
    if (sum(is.na(output$keep)) > 0) {
      output$keep[is.na(output$keep)] <- F
    }
  }
  return(output)
}

#' Choose best replicates to keep
#'
#' Reduce sample list by subsetting replicates based on a specified metric.
#'
#' Detects replicate samples from a sample list using ID duplicates and retains only replicate(s)
#' that meet selection criteria for the user-specified metric.
#'
#' @param data Data frame with at least two columns, one with a sample ID for detecting replicates
#' and another with a metric to use for selecting replicates.
#' @param id Sample ID to use for detecting replicates.
#' @param metric Metric name to use for selecting replicates.
#' @param max_replicates Maximum number of replicates to keep from a replicate set.
#' @param decreasing If 'TRUE', choose replicates by starting with the highest metric value.
#' Otherwise choose replicates by starting with the lowest metric value.
#' @return A vector of indices for selected sample replicates as well as samples without
#' replicates. Indices correspond to rows in `data`.
#' @export
select_replicates <- function(data, id, metric, max_replicates = 1, decreasing = T) {
  if (!id %in% colnames(data)) {
    stop("Error: ID column name not found in 'data'")
  }
  if (!metric %in% colnames(data)) {
    stop("Error: metric column name not found in 'data'")
  }
  if ("initial_index" %in% colnames(data)) {
    stop("Error: 'row_index' is an existing column name in 'data'. Please rename this column")
  }
  if (max_replicates < 1) {
    stop("Error: 'max_replicates' less than 1")
  }
  # Circumvent replicate detection if no replicates exist
  if (length(unlist(data[, id])) == length(unique(unlist(data[, id])))) {
    return(seq_len(nrow(data)))
  }
  data_tmp <- data.frame(row_index = seq_len(nrow(data)), data)
  id_counts <- table(unlist(data_tmp[, id]))
  replicates_list <- sapply(
    names(id_counts)[id_counts > 1],
    function(x) {
      data_tmp[which(data_tmp[, id] == x), ]
    },
    simplify = F
  )
  replicates_list_subset <- sapply(
    replicates_list,
    function(replicates_df) {
      sorted_index <- order(as.numeric(replicates_df[, metric]), decreasing = decreasing)
      replicates_df[sorted_index[1:max_replicates], ]
    },
    simplify = F
  )
  index_set1 <- do.call(rbind, replicates_list_subset)[, "row_index"]
  index_set2 <- setdiff(data_tmp[, "row_index"], do.call(rbind, replicates_list)[, "row_index"])
  selected_indices <- sort(c(index_set1, index_set2))
  return(selected_indices)
}
