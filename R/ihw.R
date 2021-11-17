# Tools for working with independent hypothesis weighting (IHW) data.
# Author: Bryan Quach (bryancquach@gmail.com)

#' Plot p-value vs. IHW variable ranking
#'
#' Plots feature p-values vs. feature IHW variable rankings.
#'
#' This plot provides insight into how informative the selected IHW variable is in relation to
#' the hypothesis being tested (e.g., differential gene expression). As IHW variable ranking
#' increases (i.e., the rank gets larger) the p-values should trend smaller for more informative
#' variables.
#'
#' @param ihw_data An ihwResult object.
#' @param truncate_y If 'TRUE', remove -log10 p-value outliers from the plot (exceed 1.5*IQR above 
#'   the third quartile).
#' @param point_color Scatterplot point fill color.
#' @param point_alpha Scatterplot point fill color alpha value.
#' @param point_size Scatterplot point size.
#' @return A ggplot object.
#' @family ihw
#' @export
ihw_pval_scatter <- function(ihw_data,
                             truncate_y = F,
                             point_color = "gray20",
                             point_alpha = 0.5,
                             point_size = 2) {
  ihw_table <- as.data.frame(ihw_data)
  if (truncate_y) {
    log_p <- -log10(ihw_table$pvalue)
    outlier_cutoff <- quantile(log_p, probs = c(0.75)) + 1.5 * IQR(log_p, na.rm = T)
    keepers <- which(log_p <= outlier_cutoff)
    ihw_table <- ihw_table[keepers, ]
  }
  output_plot <- ggplot(ihw_table, aes(x = rank(covariate, na.last = F), y = -log10(pvalue))) +
    geom_point(col = point_color, alpha = point_alpha, size = point_size) +
    labs(x = "Variable rank by increasing value", y = "-log10(p)") +
    theme(
      plot.margin = unit(c(0.5,0.5,0.5,1), units = "cm"),
      title = element_text(size = 18),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18),
      axis.title.y = element_text(vjust = 3),
      axis.title.x = element_text(vjust = -1)
    )
  return(output_plot)
}

#' Stratified p-value density distributions
#'
#' Plot p-value density distribution for IHW variable strata.
#'
#' This plot partitions the p-values by quantiles of the IHW variable to generate multiple p-value 
#' density distributions. For an informative IHW variable, the larger the quantile ID (strata) the
#' more skew there should be towards enrichment of smaller p-values.
#'
#' @param ihw_data An ihwResult object.
#' @param num_quantiles Number of strata into which the independent variable is divided.
#' @param fill_color Fill color for area under density distributions.
#' @param fill_alpha Fill color alpha value for area under density distributions.
#' @return A ggplot object.
#' @family ihw
#' @export
ihw_pval_ridges <- function(ihw_data, num_quantiles = 10, fill_color = "gray70", fill_alpha = 1) {
  ihw_table <- as.data.frame(ihw_data)
  breakpoints <- quantile(ihw_table$covariate, probs = seq(0, 1, by = 1 / num_quantiles))
  ihw_table$stratum <- cut(ihw_table$covariate, breaks = breakpoints, labels = 1:num_quantiles)
  output_plot <- ggplot(ihw_table, aes(x = -log10(pvalue), y = stratum)) +
    stat_density_ridges(
      quantile_lines = T,
      quantiles = 2,
      fill = fill_color,
      alpha = fill_alpha,
      size = 1,
      na.rm = T
    ) +
    labs(x = "-log10(p)", y = "Stratum") +
    theme_ridges(center_axis_labels = 1) +
    theme(
      legend.position = "none",
      plot.margin = unit(c(0.5, 0.25, 0.25, 0.25), "cm"),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18, hjust = 0, vjust = 0)
    )
  return(output_plot)
}

#' Stratified p-value histograms
#'
#' Plot p-value histograms for IHW variable strata.
#'
#' This plot partitions the p-values into equal sized strata of the IHW variable to generate 
#' multiple p-value histograms. For an informative IHW variable, the larger the stratum ID the
#' more skew there should be towards enrichment of smaller p-values.
#'
#' @param ihw_data An ihwResult object.
#' @param num_strata Number of strata into which the independent variable is divided.
#' @param fill_color Fill color for area under density distributions.
#' @param fill_alpha Fill color alpha value for area under density distributions.
#' @return A ggplot object.
#' @family ihw
#' @export
ihw_stratified_pval_hist <- function(ihw_data,
                                     num_strata = nlevels(groups_factor(ihw_data)),
                                     fill_color = "gray20",
                                     fill_alpha = 1) {
  ihw_table <- as.data.frame(ihw_data)
  ihw_table$stratum <- groups_by_filter(ihw_table$covariate, nbins = num_strata)
  output_plot <- ggplot(ihw_table, aes(x = pvalue)) +
    labs(title = "Stratified p-value histograms", x = "p-value", y = "Frequency") +
    geom_histogram(
      binwidth = 0.05, 
      boundary = 0, 
      color = "white", 
      fill = fill_color, 
      alpha = fill_alpha
    ) +
    facet_wrap(~stratum, nrow = ceiling(num_strata / 4)) +
    theme(
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = "cm"),
      title = element_text(size = 18),
      text = element_text(size = 18),
      axis.title = element_text(size = 18),
      axis.title.y = element_text(vjust = 3),
      axis.title.x = element_text(vjust = -2),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
      panel.grid = element_line(size = 1)
    )
  return(output_plot)
}

#' Stratified p-value eCDFs
#' 
#' Plot empirical CDFs for p-values stratified by the IHW variable.
#'
#' This plot partitions p-values into multiple eCDFs by IHW variable strata. An informative 
#' IHW variable will show a strongly increasing area under curve as the stratum ID increases.
#'
#' @param ihw_data An ihwResult object.
#' @param num_strata Number of strata into which the independent variable is divided.
#' @return A ggplot object.
#' @family ihw
#' @export
ihw_pval_cdf <- function(ihw_data,
                         num_strata = nlevels(groups_factor(ihw_data))) {
  ihw_table <- as.data.frame(ihw_data)
  ihw_table$stratum <- groups_by_filter(ihw_table$covariate, nbins = num_strata)
  output_plot <- ggplot(ihw_table, aes(x = pvalue, col = stratum)) +
    stat_ecdf(geom = "step", size = 1.5) +
    labs(color = "Stratum", x = "p-value", y = "Proportion") +
    theme(
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      title = element_text(size = 18),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.title = element_text(size = 18),
      axis.title.y = element_text(vjust = 3),
      axis.title.x = element_text(vjust = -1),
      panel.grid = element_line(size = 1)) +
      guides(colour = guide_legend(override.aes = list(size = 10)))
  return(output_plot)
}

#' Stratified IHW p-value weights
#'
#' Plot p-value weights for each strata partitioned by fold.
#'
#' Produces either a barplot or line plot depending on if the IHW variable is nominal or ordinal.
#' A stable p-value weighting scheme is expected to show consistent weights across folds within a 
#' stratum and increasing weights with increasing stratum ID.
#'
#' @param ihw_data An ihwResult object.
#' @param scale Scale of the IHW variable. Either 'nominal' or 'ordinal'.
#' @param point_size Point size for ordinal plot.
#' @param line_size Line size for ordinal plot.
#' @param alpha Plot colors alpha value.
#' @return A ggplot object.
#' @family ihw
#' @export
ihw_partitioned_pval_weights <- function(ihw_data,
                                         scale = covariate_type(ihw_data),
                                         point_size = 2,
                                         line_size = 1.5,
                                         alpha = 1) {
  fold_weights <- weights(ihw_data, levels_only = T)
  fold_ids <- factor(seq_len(nfolds(ihw_data)))
  plot_data <- expand.grid(
    stratum = seq_len(nlevels(groups_factor(ihw_data))),
    fold = fold_ids
  )
  plot_data$weight <- mapply(
    function(x, y) {
      return(fold_weights[x, y])
    },
    plot_data$stratum,
    plot_data$fold
  )
  output_plot <- switch(
    scale,
    nominal = {
      ggplot(plot_data, aes(x = fold, y = weight, fill = fold)) +
        geom_bar(stat = "identity", position = "dodge", alpha = alpha, color = "white") +
        labs(title = "P-value weighting by stratum", y = "Weight", fill = "Fold") +
        facet_wrap(~stratum) +
        theme(
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = "cm"),
          title = element_text(size = 18),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          text = element_text(size = 18),
          axis.title = element_text(size = 18),
          axis.title.y = element_text(vjust = 3),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text = element_text(size = 18),
          axis.text.x = element_blank(),
          panel.grid = element_line(size = 1)
        )
     },
     ordinal = {
       ggplot(plot_data, aes(x = stratum, y = weight, col = fold))+
        geom_point(size = point_size, alpha = alpha) +
        geom_line(size = line_size, alpha = alpha) +
        labs(x = "Stratum", y = "P-value weight", color = "Fold") +
        scale_x_continuous(breaks = scales::pretty_breaks()) +
        theme(
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = "cm"),
          title = element_text(size = 18),
          text = element_text(size = 18),
          axis.title = element_text(size = 18),
          axis.title.y = element_text(vjust = 3),
          axis.title.x = element_text(vjust = -2),
          axis.text = element_text(size = 18),
          axis.text.x = element_text(vjust = -1),
          panel.grid = element_line(size = 1)
        )
     },
     stop(sprintf("Invalid scale = '%s'", scale))
  )
  return(output_plot)
}

#' Plot IHW weight decision boundaries
#'
#' Plot p-value weight decision boundaries by strata partitioned by fold.
#'
#' @param ihw_data An ihwResult object.
#' @param scale Scale of the IHW variable. Either 'nominal' or 'ordinal'.
#' @param line_size Line size for decision boundary.
#' @param alpha Plot colors alpha value.
#' @return A ggplot object.
#' @family ihw
#' @export
ihw_pval_decision_boundary <- function(ihw_data,
                                       scale = covariate_type(ihw_data),
                                       line_size = 1.5,
                                       alpha = 1){
  if (scale != "ordinal") {
    warnings("Decision boundary plot only implemented for 'ordinal' scale")
    return(NA)
  }
  fold_weights <- weights(ihw_data, levels_only = T)      
  fold_ids <- factor(seq_len(nfolds(ihw_data)))
  plot_data <- expand.grid(
    stratum = seq_len(nlevels(groups_factor(ihw_data))),
    fold = fold_ids
  )
  plot_data$weight <- mapply(
    function(x, y) {
      return(fold_weights[x, y])
    },
    plot_data$stratum,
    plot_data$fold
  )
  ihw_table <- as.data.frame(ihw_data)
  ihw_table <- ihw_table[! is.na(ihw_table$pvalue), ]
  min_pval <- 10 ^ (-ceiling(log10(nrow(ihw_table))) - 6)
  if (sum(ihw_table$pvalue < min_pval) > 0) {
    ihw_table$pvalue[ihw_table$pvalue < min_pval] <- min_pval
  }
  ihw_table$ranking <- rank(ihw_table$covariate, ties.method = "min")
  min_rank <- by(
    ihw_table,
    ihw_table$group,
    function(df) {
      return(min(df$ranking))
    }
  )
  p_threshold <- by(
    ihw_table,
    list(ihw_table$group, ihw_table$fold),
    function(df) {
      passed <- (df$adj_pvalue <= ihw_data@alpha)
      if (! any(passed)) {
        return(min_pval)
      } else {
        return(max(df$pvalue[passed], na.rm = T))
      }
    }
  )
  plot_data$min_rank <- min_rank[as.character(plot_data$stratum)]
  threshold_indices <- cbind(as.character(plot_data$stratum), as.character(plot_data$fold))
  plot_data$p_threshold <- p_threshold[threshold_indices]
  # Add another bin for plotting
  extra_data <- subset(plot_data, stratum == nlevels(groups_factor(ihw_data)))
  extra_data$min_rank <- rep(nrow(ihw_table), nrow(extra_data))
  plot_data <- rbind(plot_data, extra_data)
  output_plot <- ggplot(ihw_table, aes(x = rank(covariate), y = -log10(pvalue))) +
    geom_hex(bins = 100, alpha = alpha) +
    geom_step(
      aes(x = min_rank, y = -log10(p_threshold), col = fold),
      data = plot_data,
      direction = "hv",
      size = line_size,
      alpha = alpha
    ) +
    labs(title = "Decision boundary", x = "Variable rank", y = "-log10(p)", color = "Fold") +
    theme(
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      title = element_text(size = 18),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18),
      axis.title.y = element_text(vjust = 3),
      axis.title.x = element_text(vjust = -1),
      panel.grid = element_line(size = 1)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 10)))
  return(output_plot)
}

#' Compare p-value transformations
#'
#' Plot IHW adjusted p-values to raw p-values or adjusted p-values using other methods.
#'
#' Creates a scatterplot that compares IHW adjusted p-values to p-values adjusted using a 
#' different method (or simply the unadjusted p-values).
#'
#' @param ihw_data An ihwResult object.
#' @param p_adjust A p-value adjustment compatible with the 'method' parameter in 'p.adjust'. If
#'   'NULL' then no p-value adjustment is made.
#' @param point_size Point size for plot.
#' @param alpha Plot colors alpha value.
#' @return A ggplot object.
#' @family ihw
#' @export
ihw_pval_comparison <- function(ihw_data,
                                p_adjust = NULL,
                                point_size = 1.5,
                                alpha = 1){
  ihw_table <- as.data.frame(ihw_data)
  if (! is.null(p_adjust)) {
    ihw_table$reference_pvalue <- p.adjust(ihw_table$pvalue, method = p_adjust)
  } else {
    ihw_table$reference_pvalue <- ihw_table$pvalue
  }
  output_plot <- ggplot(ihw_table, aes(x = reference_pvalue, y = adj_pvalue, col = group)) +
    geom_point(size = point_size) +
    #scale_colour_hue(l = 70, c = 150, drop = FALSE) +
    labs(color = "Stratum", x = "Reference p-value", y = "IHW adjusted p-value") +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      title = element_text(size = 18),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.title = element_text(size = 18),
      axis.title.y = element_text(vjust = 3),
      axis.title.x = element_text(vjust = -1),
      panel.grid = element_line(size = 1)) +
    guides(colour = guide_legend(override.aes = list(size = 5)))
  return(output_plot)
}

#' IHW diagnostic plots
#'
#' Create diagnostic plots for independent hypothesis weighting.
#'
#' @description
#' Creates ggplot objects from an ihwResult object. The ggplot objects correspond to eight plots
#' that assess if IHW is a reasonable method to use with the user-specified independent variable:
#'   * log p-value vs. independent variable ranking
#'   * log p-value density plots stratified by independent variable
#'   * p-value histograms stratified by independent variable
#'   * p-value cumulative distributions stratified by independent variable
#'   * mean IHW weight stratified by independent variable for multiple folds of the data
#'   * decision boundary relative to p-values and independent variable ranking
#'   * IHW adjusted p-values vs. raw p-values
#'   * IHW adjusted p-values vs. adjusted p-values using a different method
#' This function is a convenient wrapper to quickly generate several diagnostic plots, but each
#' plot can be generated separately using its respective function.
#'
#' @param ihw_data An ihwResult object.
#' @return A list of ggplot objects.
#' @family ihw
#' @seealso \code{\link{ihw_pval_scatter}}, \code{\link{ihw_pval_ridges}}, 
#'   \code{\link{ihw_stratified_pval_hist}}, \code{\link{ihw_pval_cdf}}, 
#'   \code{\link{ihw_partitioned_pval_weights}}, \code{\link{ihw_pval_comparison}}, 
#'   \code{\link{ihw_pval_decision_boundary}}, 
#' @export
ihw_diagnostic_plots <- function(ihw_data, num_strata = nlevels(groups_factor(ihw_data))) {
  plot_list <- list()
  ihw_table <- as.data.frame(ihw_data)
  plot_list$pval_scatter <- ihw_pval_scatter(ihw_data)
  plot_list$pval_ridges <- ihw_pval_ridges(ihw_data)
  plot_list$stratified_pval_hist <- ihw_stratified_pval_hist(ihw_data)
  plot_list$pval_cdf <- ihw_pval_cdf(ihw_data)
  plot_list$partitioned_pval_weights <- ihw_partitioned_pval_weights(ihw_data)
  plot_list$weighting_decision_boundaries <- ihw_pval_decision_boundary(ihw_data)
  plot_list$pval_comparison_raw <- ihw_pval_comparison(ihw_data, p_adjust = NULL)
  plot_list$pval_comparison_adjusted <- ihw_pval_comparison(
    ihw_data,
    p_adjust = ihw_data@adjustment_type
  )
  return(plot_list)
}


