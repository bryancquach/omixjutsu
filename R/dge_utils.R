# Generic tools for differential gene expression (DGE) analysis
# Author: Bryan Quach (bryancquach@gmail.com)

#' Merge tximport objects
#'
#' Combines multiple tximport objects into a single tximport object.
#'
#' Merges tximport objects to create a single tximport object with a common subset of genes across
#' all of the input objects. Also creates a categorical variable to distinguish the source object
#' of each sample. Requires all row and sample names to be unique.
#'
#' @param txi_list A list where each element is a tximport object.
#' @return A list with the following components:
#' * `txi` A tximport object.
#' * `sample_source` A categorical variable for each sample that corresponds to the source
#' tximport object for the sample.
#' @export
merge_txi <- function(txi_list) {
  txi_data <- list(txi = NULL, sample_source = NULL)
  if (length(txi_list) < 2) {
    warning("Less than 2 tximport objects provided in list")
    txi_data$txi <- txi_list[[1]]
    txi_data$sample_source <- rep(1, ncol(txi_list[[1]]$counts))
    names(txi_data$sample_source) <- colnames(txi_list[[1]]$counts)
    return(txi_data)
  }
  counts_from_abundance <- sapply(
    txi_list,
    function(txi) {
      return(txi$countsFromAbundance)
    }
  )
  counts_from_abundance <- unique(counts_from_abundance)
  if (length(counts_from_abundance) != 1) {
    warning("'CountsFromAbundance' values not consistent across objects")
  }
  has_dup_samples <- sapply(
    txi_list,
    function(txi) {
      has_dups <- any(
        c(
          ncol(txi$counts) != length(unique(colnames(txi$counts))),
          ncol(txi$abundance) != length(unique(colnames(txi$abundance))),
          ncol(txi$length) != length(unique(colnames(txi$length)))
        )
      )
      return(has_dups)
    },
    simplify = F
  )
  if (any(unlist(has_dup_samples))) {
    stop("Error: column names not unique within a tximport object")
  }
  has_inconsistent_samples <- sapply(
    txi_list,
    function(txi) {
      sample_ids <- c(colnames(txi$counts), colnames(txi$abundance), colnames(txi$length))
      id_freq <- table(sample_ids)
      return(any(id_freq < 3))
    },
    simplify = F
  )
  if (any(unlist(has_inconsistent_samples))) {
    stop("Error: column names not consistent within a tximport object")
  }
  txi_list <- sapply(
    txi_list,
    function(txi) {
      sorted_txi <- txi
      sorted_txi$counts <- txi$counts[order(rownames(txi$counts)), ]
      sorted_txi$abundance <- txi$abundance[order(rownames(txi$abundance)), ]
      sorted_txi$length <- txi$length[order(rownames(txi$length)), ]
      return(sorted_txi)
    },
    simplify = F
  )
  # Counts processing
  feature_names_list <- sapply(
    txi_list,
    function(txi) {
      return(rownames(txi$counts))
    }
  )
  are_unique_ids <- sapply(
    feature_names_list,
    function(ids) {
      return(length(ids) == length(unique(ids)))
    }
  )
  if (!all(unlist(are_unique_ids))) {
    stop("Error: Row names in 'counts' must all be unique")
  }
  feature_freq <- table(unlist(feature_names_list))
  if (any(feature_freq < length(txi_list))) {
    warning("Row names not consistent across 'counts' elements")
  }
  intersection_index <- which(feature_freq == length(txi_list))
  intersection_ids <- names(feature_freq)[intersection_index]
  counts_list <- sapply(
    txi_list,
    function(txi) {
      txi$counts[intersection_ids, ]
    },
    simplify = F
  )
  counts_merged <- do.call(cbind, counts_list)
  # Abundance processing
  feature_names_list <- sapply(
    txi_list,
    function(txi) {
      return(rownames(txi$abundance))
    }
  )
  are_unique_ids <- sapply(
    feature_names_list,
    function(ids) {
      return(length(ids) == length(unique(ids)))
    }
  )
  if (!all(unlist(are_unique_ids))) {
    stop("Error: Row names in 'abundance' must all be unique")
  }
  feature_freq <- table(unlist(feature_names_list))
  if (any(feature_freq < length(txi_list))) {
    warning("Row names not consistent across 'abundance' elements")
  }
  intersection_index <- which(feature_freq == length(txi_list))
  intersection_ids <- names(feature_freq)[intersection_index]
  abundance_list <- sapply(
    txi_list,
    function(txi) {
      txi$abundance[intersection_ids, ]
    },
    simplify = F
  )
  abundance_merged <- do.call(cbind, abundance_list)
  # Length processing
  feature_names_list <- sapply(
    txi_list,
    function(txi) {
      return(rownames(txi$length))
    }
  )
  are_unique_ids <- sapply(
    feature_names_list,
    function(ids) {
      return(length(ids) == length(unique(ids)))
    }
  )
  if (!all(unlist(are_unique_ids))) {
    stop("Error: Row names in 'length' must all be unique")
  }
  feature_freq <- table(unlist(feature_names_list))
  if (any(feature_freq < length(txi_list))) {
    warning("Row names not consistent across 'length' objects")
  }
  intersection_index <- which(feature_freq == length(txi_list))
  intersection_ids <- names(feature_freq)[intersection_index]
  length_list <- sapply(
    txi_list,
    function(txi) {
      txi$length[intersection_ids, ]
    },
    simplify = F
  )
  length_merged <- do.call(cbind, length_list)
  # Create labels for sample txi source
  txi_num_samples <- sapply(
    txi_list,
    function(txi) {
      return(ncol(txi$counts))
    },
    simplify = F
  )
  sample_labels <- sapply(
    seq_len(length(txi_num_samples)),
    function(label) {
      rep(x = label, times = txi_num_samples[[label]])
    }
  )
  sample_labels <- unlist(sample_labels)
  # Create new tximport object
  final_txi <- list(
    abundance = abundance_merged,
    counts = counts_merged,
    length = length_merged,
    countsFromAbundance = counts_from_abundance
  )
  txi_data$txi <- final_txi
  txi_data$sample_source <- sample_labels
  names(txi_data$sample_source) <- colnames(txi_data$txi$counts)
  return(txi_data)
}

#' Subsets the tables of a tximport object
#'
#' Subsets tximport tables using specified IDs.
#'
#' Uses a user-specified list of row or column names to subset the `counts`, `abundance`, and
#' `length` elements of the tximport object. Orders the subset according to the order in the
#' user-specified list.
#'
#' @param txi The tximport object to subset.
#' @param ids The list of row or column indices/names to use for subsetting.
#' @param dimensions The dimenions to subset. `1` for rows and `2` for columns.
#' @return A tximport object.
#' @export
subset_txi <- function(txi, ids, dimension) {
  if (!dimension %in% c(1, 2)) {
    stop("Error: 'dimension' must be '1' (row) or '2' (column)")
  }
  index <- ids
  if (is.factor(index)) {
    index <- as.character(index)
  }
  if (is.character(index)) {
    if (dimension == 1) {
      if (sum(rownames(txi$counts) %in% index) != length(index)) {
        stop("IDs missing from txi object")
      }
      index <- match(x = index, table = rownames(txi$counts))
    } else {
      if (sum(colnames(txi$counts) %in% index) != length(index)) {
        stop("IDs missing from txi object")
      }
      index <- match(x = index, table = colnames(txi$counts))
    }
  }
  txi_subset <- txi
  if (dimension == 1) {
    txi_subset$counts <- txi_subset$counts[index, ]
    txi_subset$abundance <- txi_subset$abundance[index, ]
    txi_subset$length <- txi_subset$length[index, ]
  } else {
    txi_subset$counts <- txi_subset$counts[, index]
    txi_subset$abundance <- txi_subset$abundance[, index]
    txi_subset$length <- txi_subset$length[, index]
  }
  return(txi_subset)
}

#' Generic method to remove expression features below a given threshold.
#'
#' Filters out features considered lowly expressed using a count threshold and sample proportion
#' set by the user. The sample proportion can also be estimated by calculating proportions of
#' values from a binary variable specified by the user.
#'
#' @title Filter lowly expressed features prior to DGE analysis
#' @param object A DESeqDataSet, DGEList, or matrix object from which to filter row features.
#' @param value_cutoff The minimum value required for a feature for greater than
#' `min_sample_fraction` proportion of samples to retain the feature.
#' @param min_sample_fraction The minimum proportion of samples with greater than `value_cutoff`
#' values to retain a feature.
#' @param threshold_variable Name of the binary variable to use for calculating sample fractions.
#' @param use_min_fraction If `TRUE`, uses the smaller of the sample fractions calculated from the
#' `threshold_variable`. Otherwise uses the larger of the sample fractions.
#' @param approximate If `TRUE`, round `min_sample_fraction` to a number with the hundreths digit
#' equal to '0' or '5'.
#' @export
low_expression_filter <- function(object, ...) {
  UseMethod("low_expression_filter", object)
}

#' @rdname low_expression_filter
#' @return A filtered matrix if a matrix is provided as input.
#' @export
low_expression_filter.matrix <- function(object,
                                         value_cutoff,
                                         sample_cutoff,
                                         approximate) {
  if (approximate) {
    floored_percent <- floor(sample_cutoff * 100)
    if (floored_percent %% 10 >= 5) {
      sample_cutoff <- (floor(floored_percent / 10) * 10 + 5) / 100
    } else {
      sample_cutoff <- floor(sample_cutoff * 10) / 10
    }
  }
  cat("\nUsing ", sample_cutoff, " as 'sample cutoff'\n", sep = "")
  fraction_samples_passed <- rowMeans(object > value_cutoff, na.rm = F)
  keep_index <- which(fraction_samples_passed > sample_cutoff)
  return(object[keep_index, ])
}

#' Plot log fold changes vs. log expression level.
#'
#' Plots the log2 fold change vs log10 expression level for each feature tested for
#' differential expression.
#'
#' @title Bland-Altman plot
#' @param object A DESeqResults or data.frame object.
#' @param sig_cutoff Adjusted p-value significance threshold.
#' @param sig_up_color Point color for statistically significant up-regulated features.
#' @param sig_down_color Point color for statistically significant down-regulated features.
#' @param nonsig_color Point color for features not statistical significant.
#' @param outlier_color Point color for statistically significant features flagged as outliers.
#' @param sig_alpha Point color alpha value for statistically significant features.
#' @param nonsig_alpha Point color alpha value for features not statistically significant.
#' @param sig_size Point size value for statistically significant features.
#' @param nonsig_size Point size value for features not statistically significant.
#' @param nonsig_size Point size value for outliers.
#' @return A ggplot2 object.
#' @export
ma_plot <- function(object, ...) {
  UseMethod("ma_plot", object)
}


#' @rdname ma_plot
#' @param basemean_colname A string denoting the column in `data` that contains the base mean
#' expression values.
#' @param fc_colname A string denoting the column in `data` that contains the log2 fold changes.
#' @param adj_pvalue_colname A string denoting the column in `data` that contains the adjusted
#' p-values.
#' @param outlier_index A numeric vector denoting row indices of outliers.
#' @export
ma_plot.data.frame <- function(data,
                               basemean_colname,
                               fc_colname,
                               adj_pvalue_colname,
                               outlier_index = NULL,
                               sig_cutoff = 0.05,
                               sig_up_color = "red4",
                               sig_down_color = "steelblue3",
                               nonsig_color = "gray60",
                               outlier_color = "goldenrod3",
                               sig_alpha = 0.8,
                               nonsig_alpha = 0.4,
                               sig_size = 3,
                               nonsig_size = 2,
                               outlier_size = sig_size * 2.5) {
  if (!all(outlier_index < nrow(data) & outlier_index > 0)) {
    stop("Not all outlier indices valid.")
  }
  plot_data <- data.frame(
    expression = log10(data[, basemean_colname, drop = T] + 1),
    fold_change = data[, fc_colname, drop = T],
    adj_pvalue = data[, adj_pvalue_colname, drop = T],
    sig_class = "nonsignificant",
    stringsAsFactors = F
  )
  is_upregulated <- (plot_data$adj_pvalue < sig_cutoff & plot_data$fold_change >= 0)
  is_downregulated <- (plot_data$adj_pvalue < sig_cutoff & plot_data$fold_change < 0)
  sig_outlier_index <- intersect(which(plot_data$adj_pvalue < sig_cutoff), outlier_index)
  if (any(is_upregulated)) {
    plot_data$sig_class[which(is_upregulated)] <- "up"
  }
  if (any(is_downregulated)) {
    plot_data$sig_class[which(is_downregulated)] <- "down"
  }
  if (length(sig_outlier_index) > 0) {
    plot_data$sig_class[sig_outlier_index] <- "outlier"
  }
  y_max <- max(plot_data$fold_change) * 1.1
  output_plot <- ggplot(plot_data, aes(x = expression, y = fold_change)) +
    geom_point(
      data = function(x) {
        return(x[x$sig_class == "nonsignificant", ])
      },
      shape = 16,
      col = nonsig_color,
      size = nonsig_size,
      alpha = nonsig_alpha
    ) +
    geom_point(
      data = function(x) {
        return(x[x$sig_class == "up", ])
      },
      shape = 24,
      col = sig_up_color,
      fill = sig_up_color,
      size = sig_size,
      alpha = sig_alpha
    ) +
    geom_point(
      data = function(x) {
        return(x[x$sig_class == "down", ])
      },
      shape = 25,
      col = sig_down_color,
      fill = sig_down_color,
      size = sig_size,
      alpha = sig_alpha
    ) +
    geom_point(
      data = function(x) {
        return(x[x$sig_class == "outlier", ])
      },
      shape = "*",
      col = outlier_color,
      size = outlier_size,
      alpha = sig_alpha
    ) +
    geom_hline(yintercept = 0, linetype = "solid", size = 1, color = "gray20", alpha = 0.75) +
    ylim(-1 * y_max, y_max) +
    labs(y = "log2(fold change)", x = "log10(mean expression + 1)") +
    theme(
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      title = element_text(size = 18),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18),
      axis.title.y = element_text(vjust = 3),
      axis.title.x = element_text(vjust = -1)
    )
  return(output_plot)
}

#' Plot log fold change vs log p-value.
#'
#' Plots the log2 fold change vs -log10 pvalue for features tested for differential expression.
#'
#' @title Volcano plot
#' @inheritParams ma_plot
#' @return A ggplot object.
#' @export
volcano_plot <- function(object, ...) {
  UseMethod("volcano_plot", object)
}

#' @rdname volcano_plot
#' @param fc_colname A string denoting the column in `data` that contains the log2 fold changes.
#' @param pvalue_colname A string denoting the column in `data` that contains the nominal p-values.
#' @param adj_pvalue_colname A string denoting the column in `data` that contains the adjusted
#' p-values.
#' @param outlier_index A numeric vector denoting row indices of outliers.
#' @export
volcano_plot.data.frame <- function(data,
                                    fc_colname,
                                    pvalue_colname,
                                    adj_pvalue_colname,
                                    outlier_index = NULL,
                                    sig_cutoff = 0.05,
                                    sig_up_color = "red4",
                                    sig_down_color = "steelblue3",
                                    nonsig_color = "gray60",
                                    outlier_color = "goldenrod3",
                                    sig_alpha = 0.8,
                                    nonsig_alpha = 0.4,
                                    sig_size = 3,
                                    nonsig_size = 2,
                                    outlier_size = sig_size * 2.5) {
  if (!all(outlier_index < nrow(data) & outlier_index > 0)) {
    stop("Not all outlier indices valid.")
  }
  plot_data <- data.frame(
    fold_change = data[, fc_colname, drop = T],
    log_p = -log10(data[, pvalue_colname, drop = T]),
    adj_pvalue = data[, adj_pvalue_colname, drop = T],
    sig_class = "nonsignificant",
    stringsAsFactors = F
  )
  is_upregulated <- (plot_data$adj_pvalue < sig_cutoff & plot_data$fold_change >= 0)
  is_downregulated <- (plot_data$adj_pvalue < sig_cutoff & plot_data$fold_change < 0)
  sig_outlier_index <- intersect(which(plot_data$adj_pvalue < sig_cutoff), outlier_index)
  if (any(is_upregulated)) {
    plot_data$sig_class[which(is_upregulated)] <- "up"
  }
  if (any(is_downregulated)) {
    plot_data$sig_class[which(is_downregulated)] <- "down"
  }
  if (length(sig_outlier_index) > 0) {
    plot_data$sig_class[sig_outlier_index] <- "outlier"
  }
  x_max <- max(plot_data$fold_change) * 1.1
  output_plot <- ggplot(plot_data, aes(x = fold_change, y = log_p)) +
    geom_point(
      data = function(x) {
        return(x[x$sig_class == "nonsignificant", ])
      },
      shape = 16,
      col = nonsig_color,
      size = nonsig_size,
      alpha = nonsig_alpha
    ) +
    geom_point(
      data = function(x) {
        return(x[x$sig_class == "up", ])
      },
      shape = 24,
      col = sig_up_color,
      fill = sig_up_color,
      size = sig_size,
      alpha = sig_alpha
    ) +
    geom_point(
      data = function(x) {
        return(x[x$sig_class == "down", ])
      },
      shape = 25,
      col = sig_down_color,
      fill = sig_down_color,
      size = sig_size,
      alpha = sig_alpha
    ) +
    geom_point(
      data = function(x) {
        return(x[x$sig_class == "outlier", ])
      },
      shape = "*",
      col = outlier_color,
      size = outlier_size,
      alpha = sig_alpha
    ) +
    geom_vline(xintercept = 0, linetype = "solid", size = 1, color = "gray20", alpha = 0.75) +
    xlim(-1 * x_max, x_max) +
    labs(y = "-log10(p)", x = "log2(fold change)") +
    theme(
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      title = element_text(size = 18),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18),
      axis.title.y = element_text(vjust = 3),
      axis.title.x = element_text(vjust = -1)
    )
  return(output_plot)
}

#' Plot per sample count distribution
#'
#' Plot feature count distribution summary statistics per sample.
#'
#' Plots the 25th, 50th, and 75th percentil of the feature count distribution for each sample.
#'
#' @param dds A DESeqDataSet object.
#' @param point_size Plot point size.
#' @param point_alpha Plot point alpha value.
#' @param y_lim Plot y-axis range.
#' @return A ggplot object.
#' @export
per_sample_count_distribution <- function(dds,
                                          normalized = T,
                                          point_size = 2.5,
                                          point_alpha = 1,
                                          y_lim = NULL) {
  log_counts <- log10(DESeq2::counts(dds, normalized = normalized) + 1)
  plot_data <- apply(
    log_counts,
    2,
    function(col_data) {
      quartiles <- quantile(col_data, probs = c(0.25, 0.5, 0.75))
      names(quartiles) <- c("lower_quartile", "median", "upper_quartile")
      return(quartiles)
    }
  )
  plot_data <- as.data.frame(t(plot_data))
  plot_data <- plot_data[order(plot_data$median), ]
  plot_data$x <- seq_len(nrow(plot_data))
  plot_data <- reshape2::melt(
    plot_data,
    measure.vars = c("lower_quartile", "median", "upper_quartile")
  )
  plot_data$variable <- factor(
    plot_data$variable,
    levels = c("upper_quartile", "median", "lower_quartile")
  )
  output_plot <- ggplot(
    plot_data,
    aes(x = x, y = value, shape = variable, col = variable, fill = variable)
  ) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_shape_manual(values = c(24, 21, 25)) +
    scale_fill_manual(values = c("red3", "gray20", "steelblue3")) +
    scale_color_manual(values = c("red3", "gray20", "steelblue3")) +
    labs(x = "Sample", y = "log10(counts + 1)") +
    theme(
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      legend.title = element_blank(),
      legend.text = element_text(size = 16),
      title = element_text(size = 18),
      axis.text = element_text(size = 18),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_text(size = 18),
      axis.title.y = element_text(vjust = 3),
      axis.title.x = element_text(vjust = -1)
    )
  if (!is.null(y_lim)) {
    output_plot <- output_plot + ylim(y_lim)
  }
  return(output_plot)
}
