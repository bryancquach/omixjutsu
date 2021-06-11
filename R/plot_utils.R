# Generic convenience functions for plotting and plot data munging
# Created by: Bryan Quach (bquach@rti.org)

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
                         y_title = ""){
  colnames(data) <- "value"
  ggout_list <- list()
  data_min <- min(data$value, na.rm = T)
  data_max <- max(data$value, na.rm = T)
  axis_min <- switch (
    as.character(sign(data_min)),
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
  axis_max <- switch (
    as.character(sign(data_max)),
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
      axis.title.y = element_blank())
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
#'   `colData()`. The columns in `assay` equate to data points in the PC scatterplot.
#' @param group_var A string or vector of strings for the grouping variable(s) to use from 
#'   `colData` to color points. If multiple variables are specified, they are combined into a
#'   single string to make a single new variable.
#' @param pc_x A numeric. The PC to plot on the x-axis.
#' @param pc_y A numeric. The PC to plot on the y-axis.
#' @param ntop A numeric. Specifies the top `ntop` rows ranked by decreasing variance to subset the
#'   data to prior to PCA.
#' @param center A logical. Should the data be zero-centered prior to PCA?
#' @param scale A logical. Should the data be scaled to unit variance prior to PCA?
#' @param equal_axes A logical. Should the scatterplot use the same axis limits for both axes?
#' @param point_size A numeric for the plot point size.
#' @param alpha A numeric for the plot point alpha level.
#' @param color A string for the color of plot points. Only used when `group_var` is `NULL`.
#' @param palette A string for the RColorBrewer palette name to use when `group_var` is specified.
#' @param return_data A logical. Should plot data be returned instead of a ggplot object?
#' @return A ggplot object unless `return_data` is `TRUE`, then a data frame with the
#'   user-specified PCs, grouping variable, and an attribute for the percent variace explained for
#'   each user-specified PC.
#' @seealso \code{\link{prcomp}} \code{\link{SummarizedExperiment}} \code{\link{assay}} 
#'   \code{\link{colData}}
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
                  return_data = F){
  pc_x_name <- paste0("PC", pc_x)
  pc_y_name <- paste0("PC", pc_y)
  row_vars <- rowVars(assay(data))
  keeper_rows <- order(row_vars, decreasing = T)[seq_len(min(ntop, length(row_vars)))]
  pca <- prcomp(t(assay(data)[keeper_rows, ]), center = center, scale. = scale)
  pct_var <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
  if (is.null(group_var)) {
    new_group_var <- "group"
    group_var_df <- data.frame(group = rep(1, ncol(data)))
    group <- factor(group_var_df$group)
  } else if (! all(group_var %in% names(colData(data)))) {
    stop("the argument 'group_var' should specify only column names from 'colData(data)'")
  } else {
    if (length(group_var) > 1) {
      tmp_df <- as.data.frame(colData(data)[, group_var, drop = FALSE])
      new_group_var <- "group"
      group_var_df <- data.frame(group = apply(tmp_df, 1, paste, collapse = ":"))
      group <- factor(group_var_df$group)
    } else {
      new_group_var <- group_var
      group_var_df <- as.data.frame(colData(data)[, group_var, drop = FALSE])
      group <- group_var_df[, group_var]
    }
  }
  plot_data <- data.frame(
    pca$x[, pc_x], 
    pca$x[, pc_y],
    group_var_df
  )
  rownames(plot_data) <- colnames(data)
  colnames(plot_data) <- c(pc_x_name, pc_y_name, new_group_var)
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
      plot.margin = unit(c(0.5,0.5,0.5,0.5), units = "cm"),
      title = element_text(size = 18),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18),
      axis.title.y = element_text(vjust = 3),
      axis.title.x = element_text(vjust = -1),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16)
    )
  if (equal_axes) {
    axis_min <- switch (
      as.character(sign(data_min)),
      "-1" = {
        data_min * 1.05
      },
      "1" = {
        data_min * 0.95
      },
      "0" = {
        data_min - binsize
      }
    )
    axis_max <- switch (
      as.character(sign(data_max)),
      "-1" = {
        data_max * 0.95
      },
      "1" = {
        data_max * 1.05
      },
      "0" = {
        data_max + binsize
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
#'   the categorical variable as the second column.
#' @param binsize A numeric value for the histogram bin widths.
#' @param colors A string vector. The histogram bar colors and jitter colors for each categorical
#'   variable value.
#' @param hist_alpha A numeric value for the alpha level for histogram bars.
#' @param box_fill A string. The fill color for the boxplot.
#' @param box_alpha A numeric value for the alpha level for the boxplot.
#' @param box_lwd A numeric value for the boxplot line width.
#' @param jitter_alpha A numeric for the alpha level for jittered data points.
#' @param jitter_size A numeric for the point size for jittered data points.
#' @param x_title A string denoting the x-axis title. Only added to the boxplot.
#' @param y_title A string denoting the y-axis title. Only added to the histogram.
#' @return A list of two ggplot objects `hist` and `boxplot`.
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
                          y_title = ""){
  colnames(data) <- c("value", "group")
  ggout_list <- list()
  data_min <- min(data$value, na.rm = T)
  data_max <- max(data$value, na.rm = T)
  axis_min <- switch (
    as.character(sign(data_min)),
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
  axis_max <- switch (
    as.character(sign(data_max)),
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

