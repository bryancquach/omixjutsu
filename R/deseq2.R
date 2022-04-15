# Tools for differential gene expression (DGE) analysis using DESeq2.
# Author: Bryan Quach (bryancquach@gmail.com)

#' @rdname low_expression_filter
#' @return A filtered DESeqDataSet object if a DESeqDataSet object is provided as input.
#' @export
low_expression_filter.DESeqDataSet <- function(dds,
                                               value_cutoff,
                                               min_sample_fraction = NULL,
                                               threshold_variable = NULL,
                                               use_min_fraction = T,
                                               approximate = F) {
  if (is.null(min_sample_fraction) & is.null(threshold_variable)) {
    stop("Error: either 'min_sample_fraction' or 'threshold_variable' must not be NULL")
  }
  count_matrix <- DESeq2::counts(dds, normalized = F)
  if (is.null(min_sample_fraction)) {
    pheno_data <- SummarizedExperiment::colData(dds)
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
  return(dds[rownames(filtered_count_matrix), ])
}

#' Filter lowly expressed features from a DESeq2DataSet (deprecated)
#'
#' Removes expression features below a given threshold from a DESeq2DataSet object.
#'
#' Included for backwards compatibility only. Use `low_expression_filter` instead.
#' Filters out features considered lowly expressed from a DDS object using a count threshold and
#' sample proportion set by the user. The sample proportion can also be estimated by calculating
#' proportions of values from a binary variable specified by the user.
#'
#' @param dds The DDS object from which to filter features.
#' @param count_cutoff The minimum number of counts required for a feature for greater than
#' `min_sample_fraction` proportion of samples to retain the feature.
#' @param min_sample_fraction The minimum proportion of samples with greater than `count_cutoff`
#' counts to retain a feature.
#' @param threshold_variable Name of the binary variable to use for calculating sample fractions.
#' @param use_min_fraction If `TRUE`, uses the smaller of the sample fractions calculated from the
#' `threshold_variable`. Otherwise uses the larger of the sample fractions.
#' @param approximate If `TRUE`, round `min_sample_fraction` to a whole number ending in '0' or '5'.
#' @return A DDS object.
#' @export
filter_dds_features <- function(dds,
                                count_cutoff,
                                min_sample_fraction = NULL,
                                threshold_variable = NULL,
                                use_min_fraction = T,
                                approximate = F) {
  if (is.null(min_sample_fraction) & is.null(threshold_variable)) {
    stop("Error: either 'min_sample_fraction' or 'threshold_variable' must not be NULL")
  }
  count_matrix <- DESeq2::counts(dds, normalized = F)
  if (!is.null(min_sample_fraction)) {
    sample_cutoff <- min_sample_fraction
    if (approximate) {
      floored_percent <- floor(sample_cutoff * 100)
      if (floored_percent %% 10 >= 5) {
        sample_cutoff <- (floor(floored_percent / 10) * 10 + 5) / 100
      } else {
        sample_cutoff <- floor(sample_cutoff * 10) / 10
      }
      cat("\nUsing ", sample_cutoff, " as 'min_sample_fraction'\n", sep = "")
    }
    fraction_samples_passed <- rowMeans(count_matrix > count_cutoff, na.rm = F)
    keep_index <- which(fraction_samples_passed > sample_cutoff)
    return(dds[keep_index, ])
  } else {
    pheno_data <- SummarizedExperiment::colData(dds)
    variable_freq <- table(pheno_data[, threshold_variable])
    variable_fraction <- variable_freq / nrow(pheno_data)
    if (use_min_fraction) {
      sample_cutoff <- min(variable_fraction)
    } else {
      sample_cutoff <- max(variable_fraction)
    }
    if (approximate) {
      floored_percent <- floor(sample_cutoff * 100)
      if (floored_percent %% 10 >= 5) {
        sample_cutoff <- (floor(floored_percent / 10) * 10 + 5) / 100
      } else {
        sample_cutoff <- floor(sample_cutoff * 10) / 10
      }
    }
    cat("\nUsing ", sample_cutoff, " as 'min_sample_fraction'\n", sep = "")
    fraction_samples_passed <- rowMeans(count_matrix > count_cutoff, na.rm = F)
    keep_index <- which(fraction_samples_passed > sample_cutoff)
    return(dds[keep_index, ])
  }
}

#' Get Cook's distance cutoff
#'
#' Determine Cook's distance outlier cutoff value.
#'
#' Uses differential expression regression model design to specify an F distribution and calculate
#' a value that equates to the user-specified F distribution percentile. This value can be used as
#' a Cook's distance cutoff above which features would be considered outliers.
#'
#' @param dds_fit A DESeqDataSet object after model fitting.
#' @param f_threshold Percentile of the F distribution used as the cutoff for calling outliers.
#' @return A Cook's distance value.
#' @seealso \code{\link{get_cooks_outliers}}, \code{\link{plot_cooks_per_sample}}
#' @export
get_cooks_cutoff <- function(dds_fit, f_threshold = 99) {
  if (f_threshold > 100 | f_threshold < 0) {
    stop("Error: 'f_threshold' must be between 0 and 100")
  }
  f_threshold <- f_threshold / 100
  m <- ncol(dds_fit)
  p <- ncol(model.matrix(DESeq2::design(dds_fit), SummarizedExperiment::colData(dds_fit))) - 1
  cooks_cutoff <- qf(f_threshold, p, m - p)
  return(cooks_cutoff)
}

#' Flag outliers using Cook's distance
#'
#' Flag expression features as outliers based on Cook's distance.
#'
#' Uses differential expression regression model design to specify an F distribution and calculate
#' a value that equates to the user-specified F distribution percentile. This value is used as a
#' Cook's distance cutoff above which features are considered outliers. Since a Cook's distance
#' value is calculated for every feature and sample, the max Cook's distance value across
#' samples for each feature is used to call outliers.
#'
#' @inheritParams get_cooks_cutoff
#' @return A boolean vector indicating if a feature was influenced by outlier samples.
#' @seealso \code{\link{get_cooks_cutoff}}, \code{\link{plot_cooks_per_sample}}
#' @export
get_cooks_outliers <- function(dds_fit, f_threshold = 99) {
  if (f_threshold > 100 | f_threshold < 0) {
    stop("Error: 'f_threshold' must be between 0 and 100")
  }
  max_cooks <- apply(SummarizedExperiment::assays(dds_fit)[["cooks"]], 1, max)
  cooks_cutoff <- get_cooks_cutoff(dds_fit, f_threshold)
  is_outlier <- (max_cooks > cooks_cutoff)
  names(is_outlier) <- rownames(dds_fit)
  return(is_outlier)
}

#' Plot Cook's distance distributions per sample.
#'
#' Visualizes Cook's distance distribution summary stats across all features for a sample.
#'
#' Uses Cook's distance values across all features and samples to create a plot with the 25th,
#' 50th, and 75th percentiles of the Cook's distance distrubtion for each sample.
#'
#' @param dds_fit A DESeqDataSet object after model fitting.
#' @param point_color Plot point color.
#' @param point_alpha Plot point color alpha value.
#' @param cooks_cutoff Cook's distance cutoff value for calling outliers.
#' @param y_lim y-axis lower and upper plot limit.
#' @return A ggplot object.
#' @seealso \code{\link{get_cooks_cutoff}}
#' @export
plot_cooks_per_sample <- function(dds_fit,
                                  point_size = 2.5,
                                  point_alpha = 1,
                                  cooks_cutoff = NULL,
                                  y_lim = NULL) {
  cooks_mat <- log10(SummarizedExperiment::assays(dds_fit)[["cooks"]])
  plot_data <- apply(
    cooks_mat,
    2,
    function(col_data) {
      quartiles <- quantile(col_data, probs = c(0.25, 0.5, 0.75), na.rm = T)
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
    labs(x = "Sample", y = "Cook's distance") +
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
  if (!is.null(cooks_cutoff)) {
    output_plot <- output_plot +
      geom_hline(
        yintercept = log10(cooks_cutoff),
        linetype = "dashed",
        color = "gray50",
        alpha = 0.5,
        size = 1.5
      )
  }
  if (!is.null(y_lim)) {
    output_plot <- output_plot +
      ylim(y_lim)
  }
  return(output_plot)
}

#' Compare Cook's distances to Wald statistic ranks
#'
#' Plot max Cook's distances per feature vs. Wald statistic rankings.
#'
#' Obtains the max Cook's distance for a feature across all samples and compares it to the Wald
#' statistic derived from the statistical test for differential expression. Extreme Cook's
#' distance values are truncated to the y-axis upper limit to aid with visualization.
#'
#' @param dds_fit A DESeqDataSet object after model fitting.
#' @param dds_results A DESeqResults object derived from `dds_fit`.
#' @param f_threshold Percentile of the F distribution used as the cutoff for calling outliers.
#' @param point_size Plot point size.
#' @param point_color Plot point color.
#' @param point_alpha Plot point color alpha value.
#' @return A ggplot object.
#' @export
plot_cooks_vs_wald <- function(dds_fit,
                               dds_results,
                               f_threshold = 99,
                               point_size = 3,
                               point_color = "gray30",
                               point_alpha = 0.5) {
  if (f_threshold > 100 | f_threshold < 0) {
    stop("Error: 'f_threshold' must be between 0 and 100")
  }
  f_threshold <- f_threshold / 100
  max_cooks <- apply(SummarizedExperiment::assays(dds_fit)[["cooks"]], 1, max)
  m <- ncol(dds_fit)
  p <- ncol(model.matrix(DESeq2::design(dds_fit), SummarizedExperiment::colData(dds_fit))) - 1
  cooks_cutoff <- qf(f_threshold, p, m - p)
  initial_wald_stat <- dds_results$stat
  initial_max_cooks <- apply(SummarizedExperiment::assays(dds_fit)[["cooks"]], 1, max)
  not_na <- (!is.na(initial_wald_stat))
  wald_stat <- initial_wald_stat[not_na]
  max_cooks <- initial_max_cooks[not_na]
  # Identify and truncate extreme Cook's distance values
  is_extreme <- unlist(
    sapply(
      max_cooks,
      function(x) {
        if (x > 30) {
          return(T)
        } else {
          return(F)
        }
      }
    )
  )
  truncated_max_cooks <- unlist(
    sapply(
      max_cooks,
      function(x) {
        return(min(x, 30))
      }
    )
  )
  plot_data <- data.frame(
    wald_rank = rank(wald_stat),
    max_cooks = truncated_max_cooks,
    is_extreme = is_extreme
  )
  output_plot <- ggplot(plot_data, aes(x = wald_rank, y = max_cooks)) +
    geom_point(
      mapping = aes(shape = is_extreme),
      size = point_size,
      color = point_color,
      alpha = point_alpha
    ) +
    geom_hline(
      yintercept = cooks_cutoff,
      linetype = "dashed",
      size = 1,
      color = "red3",
      alpha = 0.75
    ) +
    labs(x = "Wald statistic rank", y = "Max Cook's distance\n(Truncated axis)") +
    theme(
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      legend.position = "none",
      title = element_text(size = 18),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18),
      axis.title.y = element_text(vjust = 3),
      axis.title.x = element_text(vjust = -1)
    )
  return(output_plot)
}

#' P-value histogram for outlier features
#'
#' Plot p-value histogram for Cook's distance-based feature outliers.
#'
#' Obtains differential expression test p-values for features flagged as outliers based on a
#' Cook's distance cutoff derived from the user-specified `f_threshold`. The F distribution that
#' corresponds to `f_threshold` is parameterized by the differential exppression regression model
#' design.
#'
#' @param dds_fit A DESeqDataSet object after model fitting.
#' @param dds_results A DESeqResults object derived from `dds_fit`.
#' @param f_threshold Percentile of the F distribution used as the cutoff for calling outliers.
#' @param fill_color Plot bar color.
#' @param fill_alpha Plot bar color alpha value.
#' @return A ggplot object.
#' @export
cooks_outlier_pval_histogram <- function(dds_fit,
                                         dds_results,
                                         f_threshold = 99,
                                         fill_color = "gray30",
                                         fill_alpha = 1) {
  if (f_threshold > 100 | f_threshold < 0) {
    stop("Error: 'f_threshold' must be between 0 and 100")
  }
  f_threshold <- f_threshold / 100
  max_cooks <- apply(SummarizedExperiment::assays(dds_fit)[["cooks"]], 1, max)
  m <- ncol(dds_fit)
  p <- ncol(model.matrix(DESeq2::design(dds_fit), SummarizedExperiment::colData(dds_fit))) - 1
  cooks_cutoff <- qf(f_threshold, p, m - p)
  is_outlier <- (max_cooks > cooks_cutoff)
  plot_data <- data.frame(pvalue = NA)
  if (any(is_outlier)) {
    plot_data <- data.frame(pvalue = dds_results$pvalue[is_outlier])
  }
  output_plot <- ggplot(plot_data, aes(x = -log10(pvalue))) +
    geom_histogram(
      position = "identity",
      binwidth = 0.25,
      color = "white",
      alpha = fill_alpha,
      fill = fill_color
    ) +
    labs(title = "Cook's distance outliers", x = "-log10(p)", y = "Frequency") +
    theme(
      plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"),
      title = element_text(size = 18),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18),
      axis.title.y = element_text(vjust = 3),
      axis.title.x = element_text(vjust = -1)
    )
  return(output_plot)
}

#' @rdname ma_plot
#' @param outliers A vector of outlier IDs.
#' @export
ma_plot.DESeqResults <- function(dds_results,
                                 outliers = NULL,
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
  if (!all(outliers %in% rownames(dds_results))) {
    warnings("Not all outliers present in 'dds_results'")
  }
  dds_results <- as.data.frame(dds_results)
  outlier_index <- which(rownames(dds_results) %in% outliers)
  output_plot <- ma_plot.data.frame(
    dds_results,
    basemean_colname = "baseMean",
    fc_colname = "log2FoldChange",
    adj_pvalue_colname = "padj",
    outlier_index = outlier_index,
    sig_cutoff = sig_cutoff,
    sig_up_color = sig_up_color,
    sig_down_color = sig_down_color,
    nonsig_color = nonsig_color,
    outlier_color = outlier_color,
    sig_alpha = sig_alpha,
    nonsig_alpha = nonsig_alpha,
    sig_size = sig_size,
    nonsig_size = nonsig_size,
    outlier_size = outlier_size
  )
  return(output_plot)
}

#' @rdname volcano_plot
#' @param outliers A vector of outlier IDs.
#' @export
volcano_plot.DESeqResults <- function(dds_results,
                                      outliers = NULL,
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
  if (!all(outliers %in% rownames(dds_results))) {
    warnings("Not all outliers present in 'dds_results'")
  }
  dds_results <- as.data.frame(dds_results)
  outlier_index <- which(rownames(dds_results) %in% outliers)
  output_plot <- volcano_plot.data.frame(
    dds_results,
    fc_colname = "log2FoldChange",
    pvalue_colname = "pvalue",
    adj_pvalue_colname = "padj",
    outlier_index = outlier_index,
    sig_cutoff = sig_cutoff,
    sig_up_color = sig_up_color,
    sig_down_color = sig_down_color,
    nonsig_color = nonsig_color,
    outlier_color = outlier_color,
    sig_alpha = sig_alpha,
    nonsig_alpha = nonsig_alpha,
    sig_size = sig_size,
    nonsig_size = nonsig_size,
    outlier_size = outlier_size
  )
  return(output_plot)
}

#' @rdname per_sample_count_distribution
#' @param normalized A boolean denoting whether the plotted counts should not be normalized
#' (`FALSE`) or normalized by size/normalization factors (`TRUE`).
#' factors.
#' @export
per_sample_count_distribution.DESeqDataSet <- function(dds,
                                                       normalized = T,
                                                       point_size = 2.5,
                                                       point_alpha = 1,
                                                       y_lim = NULL) {
  log_counts <- log10(DESeq2::counts(dds, normalized = normalized) + 1)
  output_plot <- per_sample_count_distribution.matrix(
    data = log_counts,
    point_size = point_size,
    point_alpha = point_alpha,
    y_lim = y_lim
  )
  return(output_plot)
}
