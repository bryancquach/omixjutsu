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
