# Generic convenience functions for plotting and plot data munging
# Created by: Bryan Quach (bryancquach@gmail.com)

#' Histogram and jittered boxplot
#'
#' Creates a histogram and jittered boxplot.
#'
#' Creates two plots that can be organized as a single column two row multi-figure. The top row
#' plot is a histogram. The bottom row plot is a horizontal boxplot with jittered data points
#' overlayed. The x-axis of both plots are fixed to the same limits for comparability across plots.
#' To create these plots but partitioned by a categorical variable, see \code{\link{hist_boxplot2}}.
#'
#' @param data A single column data frame with numeric values to be plotted.
#' @param binsize A numeric value for the histogram bin widths.
#' @param hist_fill A string. The fill color for histogram bars.
#' @param hist_alpha A numeric value for the alpha level for histogram bars.
#' @param box_fill A string. The fill color for the boxplot.
#' @param box_alpha A numeric value for the alpha level for the boxplot.
#' @param box_lwd A numeric value for the boxplot line width.
#' @param jitter_color A string. The point color for jittered data points.
#' @param jitter_alpha A numeric for the alpha level for jittered data points.
#' @param jitter_size A numeric for the point size for jittered data points.
#' @param x_title A string denoting the x-axis title. Only added to the boxplot.
#' @param y_title A string denoting the y-axis title. Only added to the histogram.
#' @return A list of two ggplot objects `hist` and `boxplot`.
#' @seealso \code{\link{hist_boxplot2}}
#' @export
hist_boxplot <- function(data,
                         binsize = diff(range(data[, 1])) / 50,
                         hist_fill = "gray10",
                         hist_alpha = 0.75,
                         box_fill = "goldenrod",
                         box_alpha = 0.5,
                         box_lwd = 1,
                         jitter_color = "gray30",
                         jitter_alpha = 0.75,
                         jitter_size = 1.75,
                         x_title = "",
                         y_title = "") {
  colnames(data) <- "value"
  ggout_list <- list()
  data_min <- min(data$value, na.rm = T)
  data_max <- max(data$value, na.rm = T)
  axis_min <- switch(as.character(sign(data_min)),
    "-1" = {
      data_min * 1.15
    },
    "1" = {
      data_min * 0.85
    },
    "0" = {
      data_min - binsize
    }
  )
  axis_max <- switch(as.character(sign(data_max)),
    "-1" = {
      data_max * 0.85
    },
    "1" = {
      data_max * 1.15
    },
    "0" = {
      data_max + binsize
    }
  )
  ggout_list$hist <- ggplot(data, aes(x = value)) +
    geom_histogram(
      position = "identity",
      binwidth = binsize,
      alpha = hist_alpha,
      color = "white",
      fill = hist_fill
    ) +
    xlim(axis_min, axis_max) +
    labs(x = x_title, y = y_title) +
    theme(
      plot.margin = unit(c(0.5, 0.5, 0, 0.5), units = "cm"),
      plot.title = element_text(size = 18),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 18),
      axis.title.y = element_text(size = 18, vjust = 3),
    )
  ggout_list$boxplot <- ggplot(data, aes(x = value, y = 1)) +
    geom_boxplot(
      position = "identity",
      alpha = box_alpha,
      color = "black",
      fill = box_fill,
      outlier.alpha = 0,
      lwd = box_lwd
    ) +
    geom_jitter(
      shape = 16,
      color = jitter_color,
      size = jitter_size,
      alpha = jitter_alpha,
      height = 0.15
    ) +
    xlim(axis_min, axis_max) +
    labs(x = x_title, y = y_title) +
    theme(
      plot.margin = unit(c(0, 0.5, 0.5, 0.5), units = "cm"),
      plot.title = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 18),
      axis.title.x = element_text(size = 18, vjust = -1),
      axis.text.y = element_blank(),
      axis.title.y = element_blank()
    )
  return(ggout_list)
}

#' ggplot2-based PCA plot
#'
#' Create an eigenvectors scatterplot from principal component analysis.
#'
#' Conducts a principal component analysis (PCA) then produces a scatterplot using the
#' user-specified eigenvectors (ie., principal components or PCs) from the eigenvector matrix. Data
#' points can be optionally color-coded based on a user-specified variable.
#'
#' @param data A SummarizedExperiment-like object. Must be compatible with `assay()` and
#' `colData()`. The columns in `assay` equate to data points in the PC scatterplot.
#' @param group_var A string or vector of strings for the grouping variable(s) to use from
#' `colData` to color points. If multiple variables are specified, they are combined into a
#' single string to make a single new variable.
#' @param pc_x A numeric. The PC to plot on the x-axis.
#' @param pc_y A numeric. The PC to plot on the y-axis.
#' @param ntop A numeric. Specifies the top `ntop` rows ranked by decreasing variance to subset the
#' data to prior to PCA.
#' @param center A logical. Should the data be zero-centered prior to PCA?
#' @param scale A logical. Should the data be scaled to unit variance prior to PCA?
#' @param equal_axes A logical. Should the scatterplot use the same axis limits for both axes?
#' @param point_size A numeric for the plot point size.
#' @param alpha A numeric for the plot point alpha level.
#' @param color A string for the color of plot points. Only used when `group_var` is `NULL`.
#' @param palette A string for the RColorBrewer palette name to use when `group_var` is specified.
#' @param return_data A logical. Should plot data be returned instead of a ggplot object?
#' @return A ggplot object unless `return_data` is `TRUE`, then a data frame with the
#' user-specified PCs, grouping variable, and an attribute for the percent variace explained for
#' each user-specified PC.
#' @seealso \code{\link{prcomp}} \code{\link{SummarizedExperiment}} \code{\link{assay}}
#' \code{\link{colData}}
ggpca <- function(data,
                  group_var = NULL,
                  pc_x = 1,
                  pc_y = 2,
                  ntop = nrow(data),
                  center = T,
                  scale = T,
                  equal_axes = F,
                  point_size = 3,
                  alpha = 0.75,
                  color = "red4",
                  palette = "Paired",
                  return_data = F) {
  pc_x_name <- paste0("PC", pc_x)
  pc_y_name <- paste0("PC", pc_y)
  row_vars <- rowVars(SummarizedExperiment::assay(data))
  keeper_rows <- order(row_vars, decreasing = T)[seq_len(min(ntop, length(row_vars)))]
  pca <- prcomp(
    t(SummarizedExperiment::assay(data)[keeper_rows, ]),
    center = center,
    scale. = scale
  )
  pct_var <- pca$sdev^2 / sum(pca$sdev^2)
  if (is.null(group_var)) {
    group_var_df <- data.frame(group = rep(1, ncol(data)))
    group <- factor(group_var_df$group)
  } else if (!all(group_var %in% names(SummarizedExperiment::colData(data)))) {
    stop("the argument 'group_var' should specify only column names from 'colData(data)'")
  } else {
    if (length(group_var) > 1) {
      tmp_df <- as.data.frame(SummarizedExperiment::colData(data)[, group_var, drop = FALSE])
      group_var_df <- data.frame(group = apply(tmp_df, 1, paste, collapse = ":"))
      group <- factor(group_var_df$group)
    } else {
      group_var_df <- as.data.frame(SummarizedExperiment::colData(data)[, group_var, drop = FALSE])
      group <- group_var_df[, group_var]
    }
  }
  plot_data <- data.frame(
    pca$x[, pc_x],
    pca$x[, pc_y],
    group
  )
  rownames(plot_data) <- colnames(data)
  colnames(plot_data) <- c(pc_x_name, pc_y_name, "group")
  if (return_data) {
    attr(plot_data, "pct_var") <- pct_var[c(pc_x, pc_y)]
    return(plot_data)
  }
  x_title <- paste0(pc_x_name, ": ", round(pct_var[pc_x] * 100), "% variance")
  y_title <- paste0(pc_y_name, ": ", round(pct_var[pc_y] * 100), "% variance")
  data_min <- min(plot_data[, pc_x_name, drop = T], plot_data[, pc_y_name, drop = T], na.rm = T)
  data_max <- max(plot_data[, pc_x_name, drop = T], plot_data[, pc_y_name, drop = T], na.rm = T)
  ggout <- ggplot(plot_data, aes_string(x = pc_x_name, y = pc_y_name, fill = group)) +
    geom_point(size = point_size, alpha = alpha, shape = 21, color = "white") +
    labs(x = x_title, y = y_title, fill = "") +
    theme(
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = "cm"),
      title = element_text(size = 18),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18),
      axis.title.y = element_text(vjust = 3),
      axis.title.x = element_text(vjust = -1),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16)
    )
  if (equal_axes) {
    axis_min <- switch(as.character(sign(data_min)),
      "-1" = {
        data_min * 1.05
      },
      "1" = {
        data_min * 0.95
      },
      "0" = {
        data_min - (diff(c(data_min, data_max)) * 0.05)
      }
    )
    axis_max <- switch(as.character(sign(data_max)),
      "-1" = {
        data_max * 0.95
      },
      "1" = {
        data_max * 1.05
      },
      "0" = {
        data_max + (diff(c(data_min, data_max)) * 0.05)
      }
    )
    ggout <- ggout +
      xlim(axis_min, axis_max) +
      ylim(axis_min, axis_max)
  }
  if (is.null(group_var)) {
    ggout <- ggout +
      scale_fill_manual(values = color) +
      theme(legend.position = "none")
  } else {
    if (class(group) == "factor") {
      ggout <- ggout + scale_fill_brewer(palette = palette)
    } else {
      ggout <- ggout + scale_fill_distiller(palette = palette, direction = 1)
    }
  }
  return(ggout)
}

#' Histogram and jittered boxplots partitioned by categories
#'
#' Creates a histogram and jittered boxplot with color-coded, separate distributions for each value
#' of a user-specified categorical variable.
#'
#' Creates two plots that can be organized as a single column two row multi-figure. The top row
#' plot is multiple overlayed histograms from data partitioned by a variable. The bottom row plot
#' is multiple horizontal boxplots with jittered data points overlayed. The x-axis of both plots
#' are fixed to the same limits for comparability across plots.
#'
#' @param data A two-column data frame with numeric values for plotting in the first column and
#' the categorical variable as the second column.
#' @param binsize A numeric value for the histogram bin widths.
#' @param colors A string vector. The histogram bar colors and jitter colors for each categorical
#' variable value.
#' @param hist_alpha A numeric value for the alpha level for histogram bars.
#' @param box_fill A string. The fill color for the boxplot.
#' @param box_alpha A numeric value for the alpha level for the boxplot.
#' @param box_lwd A numeric value for the boxplot line width.
#' @param jitter_alpha A numeric for the alpha level for jittered data points.
#' @param jitter_size A numeric for the point size for jittered data points.
#' @param x_title A string denoting the x-axis title. Only added to the boxplot.
#' @param y_title A string denoting the y-axis title. Only added to the histogram.
#' @return A list of two ggplot objects `hist` and `boxplot`.
#' @seealso \code{\link{hist_boxplot}}
#' @export
hist_boxplot2 <- function(data,
                          binsize = diff(range(data[, 1])) / 25,
                          colors = NULL,
                          hist_alpha = 0.75,
                          box_fill = "gray50",
                          box_alpha = 0.5,
                          box_lwd = 1,
                          jitter_alpha = 0.75,
                          jitter_size = 1.75,
                          x_title = "",
                          y_title = "") {
  colnames(data) <- c("value", "group")
  ggout_list <- list()
  data_min <- min(data$value, na.rm = T)
  data_max <- max(data$value, na.rm = T)
  axis_min <- switch(as.character(sign(data_min)),
    "-1" = {
      data_min * 1.1
    },
    "1" = {
      data_min * 0.9
    },
    "0" = {
      data_min - binsize
    }
  )
  axis_max <- switch(as.character(sign(data_max)),
    "-1" = {
      data_max * 0.9
    },
    "1" = {
      data_max * 1.1
    },
    "0" = {
      data_max + binsize
    }
  )
  ggout_list$hist <- ggplot(data, aes(x = value, fill = group)) +
    geom_histogram(
      position = "identity",
      binwidth = binsize,
      alpha = hist_alpha,
      color = "white",
    ) +
    xlim(axis_min, axis_max) +
    labs(x = x_title, y = y_title, fill = "") +
    theme(
      plot.margin = unit(c(0.5, 0.5, 0, 0.5), units = "cm"),
      plot.title = element_text(size = 18),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 18),
      axis.title.y = element_text(size = 18, vjust = 3),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16)
    )
  ggout_list$boxplot <- ggplot(data, aes(x = value, y = group, color = group)) +
    geom_boxplot(
      position = "identity",
      alpha = box_alpha,
      color = "black",
      fill = box_fill,
      outlier.alpha = 0,
      lwd = box_lwd
    ) +
    geom_jitter(
      shape = 16,
      size = jitter_size,
      alpha = jitter_alpha,
      height = 0.15
    ) +
    xlim(axis_min, axis_max) +
    labs(x = x_title, y = y_title, color = "") +
    theme(
      plot.margin = unit(c(0, 0.5, 0.5, 0.5), units = "cm"),
      plot.title = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 18),
      axis.title.x = element_text(size = 18, vjust = -1),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16)
    )
  if (is.null(colors)) {
    ggout_list$hist <- ggout_list$hist + scale_fill_brewer(palette = "Dark2")
    ggout_list$boxplot <- ggout_list$boxplot + scale_color_brewer(palette = "Dark2")
  } else {
    ggout_list$hist <- ggout_list$hist + scale_fill_manual(values = colors)
    ggout_list$boxplot <- ggout_list$boxplot + scale_color_manual(values = colors)
  }
  return(ggout_list)
}

#' Plot p-value histogram
#'
#' Plot p-value histogram.
#'
#' @param pvalues A vector of pvalues.
#' @param bin_width Size of each histogram bin.
#' @param bin_fill Bin fill color.
#' @param alpha Bin fill color alpha value.
#' @return A ggplot object.
#' @export
pval_histogram <- function(pvalues, bin_width = 0.025, bin_fill = "gray10", alpha = 0.8) {
  plot_data <- data.frame(pvalue = pvalues)
  output_plot <- ggplot(plot_data, aes(x = pvalue)) +
    geom_histogram(
      position = "identity",
      binwidth = bin_width,
      alpha = alpha,
      color = "white",
      fill = bin_fill
    ) +
    xlim(0, 1) +
    labs(y = "Frequency", x = "p-value") +
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

#' Quantile-quantile plot
#'
#' Plots observed vs. expected p-values.
#'
#' Plot observed p-values vs. expected p-values. Expected p-values are assumed to follow a
#' uniform distribution.
#'
#' @param pvalues A vector of pvalues with associated feature IDs.
#' @param outliers A vector of outlier IDs that correspond to the names in 'pvalues'.
#' @param sig_cutoff Adjusted p-value significance threshold.
#' @param plot_lambda If 'TRUE' calculate the genomic inflation factor and overlay it on the plot.
#' @param df Degrees of freedom on the theoretical distribution. Used in calculating the genomic
#' inflation factor. Only relevant when `plot_lambda` is 'TRUE'.
#' @return A ggplot object.
#' @export
pval_qqplot <- function(pvalues, outliers = NULL, sig_cutoff = 0.05, plot_lambda = T, df = 1) {
  if (!all(outliers %in% names(pvalues))) {
    warnings("Not all outliers present in 'pvalues'")
  }
  plot_data <- data.frame(
    pvalues = pvalues,
    log_p = -log10(pvalues),
    is_outlier = (names(pvalues) %in% outliers)
  )
  sorted_p <- sort(plot_data$pvalues[which(!plot_data$is_outlier)])
  if (plot_lambda) {
    lambda <- qchisq(median(sorted_p), df, lower.tail = F) / qchisq(0.5, df, lower.tail = F)
    lambda <- round(lambda, 3)
  }
  keepers <- which(plot_data$log_p < Inf & (!plot_data$is_outlier))
  plot_data <- plot_data[keepers, , drop = F]
  plot_data <- plot_data[order(plot_data$log_p, decreasing = T), ]
  num_pvals <- nrow(plot_data)
  plot_data$expected_log_p <- -log10((1:num_pvals) / (num_pvals + 1))
  plot_data$ci_lower <- -log10(
    qbeta(
      sig_cutoff / 2,
      1:num_pvals,
      (num_pvals + 1) - (1:num_pvals)
    )
  )
  plot_data$ci_upper <- -log10(
    qbeta(
      1 - (sig_cutoff / 2),
      1:num_pvals,
      (num_pvals + 1) - (1:num_pvals)
    )
  )
  plot_data$ci_expected <- -log10(((1:num_pvals) - 0.5) / num_pvals)
  xy_max <- max(plot_data$expected_log_p, plot_data$log_p) + 1
  output_plot <- ggplot(plot_data, aes(x = expected_log_p, y = log_p)) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "solid",
      size = 1.5,
      color = "red3",
      alpha = 0.75
    ) +
    geom_point(size = 2, shape = 16, color = "gray30", alpha = 0.75) +
    geom_line(
      aes(x = ci_expected, y = ci_lower),
      size = 1.5,
      linetype = "dashed",
      color = "gray60",
      alpha = 0.75
    ) +
    geom_line(
      aes(x = ci_expected, y = ci_upper),
      size = 1.5,
      linetype = "dashed",
      color = "gray60",
      alpha = 0.75
    ) +
    xlim(0, xy_max) +
    ylim(0, xy_max) +
    labs(
      x = expression(-log ~ ""["10"] ~ "(Expected p-value)"),
      y = expression(-log ~ ""["10"] ~ "(Observed p-value)")
    ) +
    theme(
      plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"),
      legend.position = "none",
      title = element_text(size = 18),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18),
      axis.title.y = element_text(vjust = 3),
      axis.title.x = element_text(vjust = -1)
    )
  if (plot_lambda) {
    output_plot <- output_plot +
      annotate(
        "text",
        x = xy_max - 1,
        y = 0,
        label = bquote(paste(lambda, " = ", .(lambda))),
        size = 6
      )
  }
  return(output_plot)
}

#' Grouped barplot
#'
#' Creates a grouped barplot.
#'
#' Creates a barplot for frequencies of two factor variables in which one variable is used to
#' stratify the other.
#'
#' @param data A data frame with columns from which to retrieve variables to plot.
#' @param x_var A string denoting the column name of the primary barplot variable.
#' @param y_var A string denoting the column name of the secondary barplot variable.
#' @param rm_ids A string vector of the row names to exclude prior to plotting.
#' @param fill A string vector denoting the fill colors for barplot bars.
#' @param alpha A numeric value for the alpha level for barplot bars.
#' @param lwd A numeric value for the barplot bar line width.
#' @param fill A string denoting the line color for barplot bars.
#' @param x_title A string denoting the x-axis title.
#' @param y_title A string denoting the y-axis title.
#' @param legend_title A string denoting the legend title. Should correspond to `y_var`.
#' @return A ggplot object.
#' @export
grouped_barplot <- function(data,
                            x_var,
                            y_var,
                            rm_ids = NULL,
                            fill = NULL,
                            alpha = 1,
                            lwd = 1,
                            line_color = "black",
                            x_title = NULL,
                            y_title = NULL,
                            legend_title = NULL) {
  if (!is.null(rm_ids)) {
    if (!all(rm_ids %in% rownames(data))) {
      warning("Not all elements in 'rm_ids' found in 'data'")
    }
    keepers <- which(!rownames(data) %in% rm_ids)
    data <- data[keepers, ]
  }
  plot_data <- reshape2::melt(
    table(data[, x_var], data[, y_var]),
    varnames = c("var1", "var2")
  )
  plot_data$var1 <- as.factor(plot_data$var1)
  plot_data$var2 <- as.factor(plot_data$var2)
  if (is.null(fill)) {
    num_factors <- nlevels(plot_data$var2)
    fill <- brewer.pal(n = max(3, num_factors), name = "Spectral")[1:num_factors]
  }
  output_plot <- ggplot(plot_data, aes(x = value, y = var1, fill = var2)) +
    geom_bar(
      stat = "identity",
      position = position_dodge(),
      alpha = alpha,
      lwd = lwd,
      color = line_color
    ) +
    scale_fill_manual(values = fill) +
    labs(x = x_title, y = y_title, fill = legend_title) +
    theme(
      plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"),
      title = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.key.size = unit(1, "cm"),
      legend.text = element_text(size = 16),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18),
      axis.title.y = element_text(vjust = 3),
      axis.title.x = element_text(vjust = -1)
    )
  return(output_plot)
}

#' Heatmap grid
#'
#' Creates a heatmap based on a matrix of values.
#'
#' @param data A data frame or matrix of numeric values.
#' @param digits An integer for how many digits to which to round cell values.
#' @param text_size A numeric value for the cell value text size.
#' @param legend_height A numeric for the height of the legend key in millimeters.
#' @param row_ids An optional string vector of row names to retain for plotting.
#' @param col_ids An optional string vector of column names to retain for plotting.
#' @param ggfill An object returned by the family of `scale_fill_*` functions for continuous values
#' that defines the color fill for the cells of the heatmap.
#' @param reorder_matrix A boolean indicating if the columns and rows should be reordered such that
#' the diagonal cells equate to the same x and y axis labels.
#' @param title A string for the plot title.
#' @param x_title A string for the x-axis title.
#' @param y_title A string for the y-axis title.
#' @param xlab_pos The x-axis label position. One of "top", "bottom".
#' @param ylab_pos The y-axis label position. One of "left", "right".
#' @return A ggplot2 object
#' @seealso \code{\link{assoc_matrix}}
#' @export
matrix_heatmap <- function(data,
                           digits = 2,
                           text_size = 4.5,
                           legend_height = 20,
                           row_ids = NULL,
                           col_ids = NULL,
                           ggfill = scale_fill_gradient2(
                             name = "",
                             low = "steelblue4",
                             mid = "white",
                             high = "red4",
                             breaks = seq(-1, 1, 0.1),
                             limits = c(-1, 1)
                           ),
                           reorder_matrix = T,
                           title = NULL,
                           x_title = NULL,
                           y_title = NULL,
                           xlab_pos = c("bottom", "top"),
                           ylab_pos = c("left", "right")) {
  if (!is.null(row_ids)) {
    if (!all(row_ids %in% rownames(data))) {
      stop("Error: Not all `row_ids` in `data`")
    } else {
      data <- data[row_ids, , drop = F]
    }
  }
  if (!is.null(col_ids)) {
    if (!all(col_ids %in% colnames(data))) {
      stop("Error: Not all `col_ids` in `data`")
    } else {
      data <- data[, col_ids, drop = F]
    }
  }
  xlab_pos <- match.arg(xlab_pos)
  ylab_pos <- match.arg(ylab_pos)
  data <- as.data.frame(data)
  measure_varnames <- colnames(data)
  data$covariate1 <- rownames(data)
  plot_data <- reshape2::melt(
    data,
    id.vars = "covariate1",
    measure.vars = measure_varnames,
    variable.name = "covariate2"
  )
  if (reorder_matrix) {
    plot_data$covariate1 <- reorder(plot_data$covariate1, -plot_data$value)
    plot_data$covariate2 <- reorder(plot_data$covariate2, plot_data$value)
  }
  output_plot <- ggplot(plot_data, aes(y = covariate1, x = covariate2, fill = value)) +
    geom_tile(colour = "white", height = 0.95, width = 0.95) +
    geom_text(aes(label = round(value, digits)), size = text_size) +
    ggfill +
    labs(x = x_title, y = y_title, title = title) +
    scale_x_discrete(position = xlab_pos) +
    scale_y_discrete(position = ylab_pos) +
    theme(
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = "cm"),
      axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 18),
      title = element_text(size = 16),
      legend.key.height = unit(legend_height, "mm"),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16)
    )
  if (xlab_pos == "top") {
    output_plot <- output_plot + theme(axis.text.x = element_text(hjust = 0))
  }
  return(output_plot)
}
