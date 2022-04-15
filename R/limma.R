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

#' Bland-Altman plot
#'
#' Plot log fold changes vs. log expression level.
#'
#' Plots the log2 fold change vs log10 expression level for each feature tested for
#' differential expression.
#'
#' @param dds_result A DESeqResults object.
#' @param outliers A vector of outlier IDs.
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
#' @return A ggplot object.
#' @export
ma_plot <- function(dds_results,
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
  plot_data <- data.frame(
    expression = log10(dds_results$baseMean),
    fold_change = dds_results$log2FoldChange,
    adj_pvalue = dds_results$padj,
    sig_class = "nonsignificant",
    stringsAsFactors = F
  )
  rownames(plot_data) <- rownames(dds_results)
  is_upregulated <- (plot_data$adj_pvalue < sig_cutoff & plot_data$fold_change >= 0)
  is_downregulated <- (plot_data$adj_pvalue < sig_cutoff & plot_data$fold_change < 0)
  is_outlier <- (plot_data$adj_pvalue < sig_cutoff & (rownames(plot_data) %in% outliers))
  if (any(is_upregulated)) {
    plot_data$sig_class[which(is_upregulated)] <- "up"
  }
  if (any(is_downregulated)) {
    plot_data$sig_class[which(is_downregulated)] <- "down"
  }
  if (any(is_outlier)) {
    plot_data$sig_class[which(is_outlier)] <- "outlier"
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
    labs(y = "log2(fold change)", x = "log10(mean expression)") +
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

#' Volcano plot
#'
#' Plot log fold change vs log p-value.
#'
#' Plots the log2 fold change vs -log10 pvalue for features tested for differential expression.
#'
#' @inheritParams ma_plot
#' @return A ggplot object.
#' @export
volcano_plot <- function(dds_results,
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
  plot_data <- data.frame(
    fold_change = dds_results$log2FoldChange,
    log_p = -log10(dds_results$pvalue),
    adj_pvalue = dds_results$padj,
    sig_class = "nonsignificant",
    stringsAsFactors = F
  )
  rownames(plot_data) <- rownames(dds_results)
  is_upregulated <- (plot_data$adj_pvalue < sig_cutoff & plot_data$fold_change >= 0)
  is_downregulated <- (plot_data$adj_pvalue < sig_cutoff & plot_data$fold_change < 0)
  is_outlier <- (plot_data$adj_pvalue < sig_cutoff & (rownames(plot_data) %in% outliers))
  if (any(is_upregulated)) {
    plot_data$sig_class[which(is_upregulated)] <- "up"
  }
  if (any(is_downregulated)) {
    plot_data$sig_class[which(is_downregulated)] <- "down"
  }
  if (any(is_outlier)) {
    plot_data$sig_class[which(is_outlier)] <- "outlier"
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
