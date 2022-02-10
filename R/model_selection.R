# Tools for model selection prior to regression modeling.
# Author: Bryan Quach (bryancquach@gmail.com)

#' Create design matrix
#'
#' Format data table to produce a design matrix.
#'
#' Creates a design matrix and prunes variables with no variance. Removes pairwise collinear
#' variables.
#'
#' @param data A data frame of potential explanatory variables for the design matrix. Rows should
#' be records and columns should be variables.
#' @param var_names A vector of column names from `data` to consider for the design matrix.
#' @param corr_cutoff A numeric value in [0,1] denoting the absolute correlation cutoff to
#' classify pairs of collinear variables.
#' @param corr_method A string indicating the type of correlation to use for determining
#' collinearity. Must be `pearson` or `spearman`.
#' @param check_binary_vars A boolean indiciating if dichotomous variables should be included in
#' the collinearity pruning.
#' @param strings_as_factors A boolean indiciating if character variables should be converted to
#' factors in the returned data frame.
#' @return A data frame of explanatory variables for regression.
#' @seealso \code{\link{get_residual_matrix}}
#' @export
get_design_mat <- function(data,
                           var_names = NULL,
                           corr_cutoff = 0.75,
                           corr_method = c("spearman", "pearson"),
                           check_binary_vars = T,
                           strings_as_factors = T) {
  if (corr_cutoff < 0 | corr_cutoff > 1) {
    stop("Error: 'corr_cutoff' must be in the range [0,1]")
  }
  corr_method <- match.arg(corr_method)
  if (length(var_names) > 0) {
    if (any(!var_names %in% colnames(data))) {
      stop("Error: Not all `var_names` found in `data`")
    }
    if (length(var_names) == 1) {
      data_subset <- data[, var_names, drop = F]
      if (length(unique(na.omit(data_subset[, 1, drop = T]))) > 1) {
        return(data_subset)
      } else {
        stop("Error: no `var_names` variable with non-zero variance")
      }
    } else {
      data <- data[, var_names]
    }
  }
  is_keeper <- apply(
    data,
    2,
    function(x) {
      num_items <- length(unique(na.omit(x)))
      return(num_items > 1)
    }
  )
  if (sum(is_keeper) > 0) {
    print(paste("Variables with non-zero variance:", sum(is_keeper)))
    if (sum(is_keeper) == 1) {
      return(data[, is_keeper, drop = F])
    }
    data <- data[, is_keeper]
    corrcheck_data <- data
  } else {
    stop("Error: no variables with non-zero variance")
  }
  # Only binary or numeric variables can be compared for collinearity
  if (check_binary_vars) {
    is_binary <- apply(
      corrcheck_data,
      2,
      function(x) {
        num_items <- length(unique(na.omit(x)))
        return(num_items == 2)
      }
    )
    print(paste("Binary variables:", sum(is_binary)))
    # Binary variables must be numeric for correlation calculations
    for (i in seq_len(length(is_binary))) {
      if (is_binary[i]) {
        corrcheck_data[, i] <- as.numeric(as.factor(corrcheck_data[, i, drop = T]))
      }
    }
  } else {
    is_binary <- rep(FALSE, ncol(corrcheck_data))
  }
  is_numeric <- sapply(corrcheck_data, class) %in% c("integer", "numeric")
  corr_candidates <- which(is_numeric | is_binary)
  print(paste("Collinearity assessment candidates:", length(corr_candidates)))
  covar_mat <- Hmisc::rcorr(as.matrix(corrcheck_data[, corr_candidates]), type = corr_method)$r
  rm_variables <- caret::findCorrelation(covar_mat, cutoff = corr_cutoff, names = T)
  print("Variables to remove:")
  print(rm_variables)
  data <- data[, !colnames(data) %in% rm_variables]
  if (strings_as_factors) {
    is_character <- sapply(data, class) == "character"
    for (i in which(is_character)) {
      data[, i] <- as.factor(data[, i, drop = T])
    }
  }
  return(data)
}

#' Calculate residual matrix
#'
#' Get residuals from fitting a linear model with a matrix of outcomes.
#'
#' Calculates regression residuals using a fixed set of explanatory variables and a matrix of
#' outcomes.
#'
#' @param design A data frame with explanatory variables.
#' @param outcomes A data frame like object with outcome variables as rows and observations as
#' columns.
#' @return A list with the following items:
#' * residuals: A numeric matrix with each column representing the residuals for an outcome. The
#' ordering follows the row order of `outcomes`.
#' * dof: Model degrees of freedom.
#' @seealso \code{\link{get_matrix_r2}} \code{\link{get_design_mat}}
#' @export
get_residual_matrix <- function(design, outcomes) {
  results <- list()
  y <- outcomes
  x <- design
  y <- t(y)
  resid_mat <- matrix(ncol = ncol(y), nrow = nrow(y))
  for (i in seq_len(ncol(y))) {
    lm_fit <- RcppEigen::fastLm(y[, i] ~ 1 + ., data = x)
    resid_mat[, i] <- lm_fit$residuals
  }
  results[["residuals"]] <- resid_mat
  results[["dof"]] <- lm_fit$df.residual
  return(results)
}

#' Calculate (adjusted) r-squared for a matrix
#'
#' Calculates regression model r-squared across a matrix of outcomes.
#'
#' @param outcomes A data frame like object with outcome variables as rows and observations as
#' columns.
#' @param residuals A numeric matrix with each column representing the residuals for an outcome.
#' The ordering follows the row order of `outcomes`.
#' @param adjusted A boolean indicating if the adjusted r-squared should be calculated in addition
#' to the unadjusted r-squared.
#' @param dof An integer denoting the degrees of freedom of the regression model used to compute
#' the residuals. Only required if computing adjusted r-squared.
#' @return A list with the following items:
#' * r2: A numeric value corresponding to the unadjusted r-squared.
#' * adjusted_r2: A numeric value (or NA) corresponding to the adjusted r-squared.
#' @seealso \code{\link{get_residual_matrix}}
#' @export
get_matrix_r2 <- function(outcomes,
                          residuals,
                          adjusted = F,
                          dof = NULL) {
  if (adjusted & is.null(dof)) {
    stop("Error: degrees of freedom required for adjusted r-squared calculation")
  }
  tss <- norm(outcomes - rowMeans(outcomes, na.rm = T), type = "F")^2
  rss <- norm(residuals, type = "F")^2
  n <- ncol(outcomes)
  r2 <- (tss - rss) / tss
  if (adjusted) {
    adjusted_r2 <- 1 - ((1 - r2) * ((n - 1) / dof))
  } else {
    adjusted_r2 <- NA
  }
  return(list(r2 = r2, adjusted_r2 = adjusted_r2))
}

#' Plot (adjusted) r-squared for multiple models
#'
#' Plots regression model r-squared vs. model size (number of variables).
#'
#' Creates a scatterplot showing the tradeoff between variance explained of the outcome
#' (e.g., gene expression) and model complexity. This function is designed to be a visualization
#' helper function for the output from `get_matrix_r2`. It allows for the plotting of only r2, only
#' adjusted r2, or both.
#'
#' @param pve_list A list of lists containing regression model information. Each list must have the
#' same element ordering, and the following elements and corresponding names are required:
#' * `r2`: The coefficient of determination. Can be calculated using `get_matrix_r2()`. Value can
#' be `NA` if `adjusted_r2` is numeric.
#' * `adjusted_r2`: The adjusted r2. Can be calculated using `get_matrix_r2()`. Value can be `NA`
#' if `r2` is numeric.
#' * `num_model_vars`: The number of variables in the model used to calculate `r2` and
#' `adjusted_r2`.
#' @param model_names A vector of model names that follows the same ordering as the lists in
#' `pve_list`. If not provided, a sequential numbering will be used.
#' @param line_color A string for the line color between points.
#' @param line_alpha A numeric for the line alpha value.
#' @param line_type A string for the line type between points.
#' @param line_size A numeric for the line size between points.
#' @param point_color A vector of colors corresponding to the groups of points (r2 vs adjusted r2).
#' If not specified, default colors from a color palette will be used.
#' @param palette A string corresponding to a color palette from which point colors are chosen.
#' @param point_alpha A numeric for the point alpha value.
#' @param point_shape A numeric for the point shape.
#' @param point_size A numeric for the point size.
#' @param label_size A numeric for the label text size..
#' @param label_line_size A numeric for the line thickness of the guide line from the label text
#' to its corresponding point.
#' @param label_line_alpha. A numeric for the line alpha value of the guide line from the label
#' text to its corresponding point.
#' @param title Plot title.
#' @param x_title X-axis title.
#' @param y_title Y-axis title.
#' @param legend_title Legend title.
#' @param xlim_multiplier A numeric constant used to increase or decrease the x-axis range.
#' @param ylim_multiplier A numeric constant used to icnrease or decrease the y-axis range.
#' @return A ggplot object.
#' @seealso \code{\link{get_matrix_r2}} \code{\link{get_residual_matrix}}
#' \code{\link{get_design_mat}}
#' @export
plot_r2 <- function(pve_list,
                    model_names = NULL,
                    line_color = "gray50",
                    line_alpha = 0.5,
                    line_type = "dashed",
                    line_size = 1.5,
                    point_color = NULL,
                    palette = "Paired",
                    point_alpha = 1,
                    point_shape = 16,
                    point_size = 5,
                    label_size = 5,
                    label_line_size = 0.75,
                    label_line_alpha = 0.5,
                    title = NULL,
                    x_title = "Model complexity (no. of variables)",
                    y_title = "Variance explained",
                    legend_title = NULL,
                    legend_position = "right",
                    xlim_multiplier = 0.25,
                    ylim_multiplier = 0.25) {
  pve_df <- data.frame(t(sapply(pve_list, unlist)))
  if (is.null(model_names)) {
    model_names <- paste0("Model ", seq_len(nrow(pve_df)))
  } else {
    if (nrow(pve_df) != length(model_names)) {
      stop("Error: `model_names` and `pve_list` are not equal length")
    }
  }
  pve_df$model <- model_names
  plot_data <- reshape2::melt(
    pve_df,
    id.vars = c("model", "num_model_vars"),
    measure.vars = c("r2", "adjusted_r2")
  )
  plot_data$label <- paste0(plot_data$model, "\n(", round(plot_data$value, 2), ")")
  # Axis limits need wiggle room for label text to fit
  x_limits <- range(plot_data$num_model_vars)
  x_buffer <- c(-1 * abs(diff(x_limits)) * xlim_multiplier, abs(diff(x_limits)) * xlim_multiplier)
  x_limits <- x_limits + x_buffer
  y_limits <- range(plot_data$value)
  y_buffer <- c(-1 * abs(diff(y_limits)) * ylim_multiplier, abs(diff(y_limits)) * ylim_multiplier)
  y_limits <- y_limits + y_buffer
  output_plot <- ggplot(
    data = plot_data,
    mapping = aes(x = num_model_vars, y = value, group = variable, label = label)
  ) +
    geom_line(linetype = line_type, color = line_color, alpha = line_alpha, size = line_size) +
    geom_point(
      aes(color = variable),
      size = point_size,
      alpha = point_alpha,
      shape = point_shape
    ) +
    xlim(x_limits[1], x_limits[2]) +
    ylim(y_limits[1], y_limits[2]) +
    ggrepel::geom_text_repel(
      box.padding = 0.75,
      max.overlaps = Inf,
      size = label_size,
      segment.alpha = label_line_alpha,
      segment.size = label_line_size
    ) +
    labs(title = title, x = x_title, y = y_title, color = legend_title) +
    theme(
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = "cm"),
      title = element_text(size = 16),
      legend.title = element_text(size = 18),
      legend.key.size = unit(1, "cm"),
      legend.text = element_text(size = 16),
      legend.text.align = 0,
      legend.position = legend_position,
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18),
      axis.title.y = element_text(vjust = 3),
      axis.title.x = element_text(vjust = -1),
    )
  if (is.null(point_color)) {
    output_plot <- output_plot +
      scale_color_brewer(palette = palette, labels = c(expression(r^2), expression("adjusted r"^2)))
  } else {
    output_plot <- output_plot +
      scale_color_manual(
        values = point_color, 
        labels = c(expression(r^2), expression("adjusted r"^2))
      )
  }
  return(output_plot)
}

#' Association matrix
#'
#' Creates an association matrix based on pairwise measures of association between variables.
#'
#' Calculates pairwise associations for a set of variables. Depending on the measure
#' of association specified, variables will be excluded if the variable type
#' (e.g., nominal or continuous) does not make sense to include in the calculations. Cramer's V
#' will only be calculated between two nominal variables. Eta-squared will only be applied to
#' nominal-continuous variable pairs. Pearson, spearman, and kendall correlations exclude nominal
#' variables with >2 values. Variable exclusions are based on the variable type as defined in
#' `data`, so these should be verified. Categorical variables with values coded as integers can be
#' mistakenly treated as continuous variables. Nominal variables should be of class `factor` (not 
#' `character`). Numeric variables should be of class `integer` or `numeric`.
#'
#' @param data A data frame with columns from which to retrieve variables to compute associations.
#' @param var_names A vector of variables names from the columns of `data` to consider.
#' @param factor_vars A vector that includes the names of variables that should be converted to 
#' factors. Must be in `data` or `var_names` if specified.
#' @param method The type of association to calculate. One of `pearson` (default), `spearman`,
#' `kendall`, `eta_squared`, `cramers_v`.
#' @param use A string giving a method for computing covariances in the presence of missing
#' values. This must be one of the strings `everything`, `all.obs`, `complete.obs`,
#' `na.or.complete`, or `pairwise.complete.obs`. See the `cor` function documentation for details
#' on what each value specifies. Only relevant for correlation metrics.
#' @return A data frame with association metric values.
#' @export
assoc_matrix <- function(data, 
                         var_names = NULL,
                         factor_vars = NULL,
                         method = c("pearson", "spearman", "kendall", "eta_squared", "cramers_v"),
                         use = c(
                           "pairwise.complete.obs", 
                           "everything", 
                           "all.obs", 
                           "complete.obs", 
                           "na.or.complete")) {
  method <- match.arg(method)
  use <- match.arg(use)
  if (length(var_names) > 0) {
    if (any(!var_names %in% colnames(data))) {
      stop("Error: Not all `var_names` found in `data`")
    }
    if (length(var_names) == 1) {
      stop("Error: `var_names` must include >=2 variables")
    }
    data <- data[, var_names]
  }
  if (!is.null(factor_vars)) {
    if(!all(factor_vars %in% colnames(data))) {
      stop("Error: Not all `factor_vars` found in `data` (or `var_names` if not NULL)")
    }
    for (i in factor_vars) {
      data[, i] <- as.factor(data[, i])
    }
  }
  print(paste("Initial variables:", ncol(data)))
  # Associations are only for variables with non-zero variance
  is_keeper <- apply(
    data,
    2,
    function(x) {
      num_items <- length(unique(na.omit(x)))
      return(num_items > 1)
    }
  )
  if (sum(is_keeper) > 1) {
    print(paste("Variables with non-zero variance:", sum(is_keeper)))
    data <- data[, is_keeper]
  } else {
    stop("Error: <2 variables with non-zero variance")
  }
  # Calculations on binary or numeric variables only
  if (method %in% c("pearson", "spearman", "kendall")) {
    is_binary <- apply(
      data,
      2,
      function(x) {
        num_items <- length(unique(na.omit(x)))
        return(num_items == 2)
      }
    )
    # Binary variables must be numeric for correlation calculations
    for (i in seq_len(length(is_binary))) {
      if (is_binary[i]) {
        data[, i] <- as.numeric(as.factor(data[, i, drop = T]))
      }
    }
    is_numeric <- sapply(data, class) %in% c("integer", "numeric")
    candidates <- which(is_numeric | is_binary)
    print(paste("Final variable set size:", length(candidates)))
    if (length(candidates) < 2) {
      stop("Error: <2 variables with non-zero variance")
    }
    data <- as.matrix(data[, candidates])
    assoc_out <- cor(data, method = method, use = use)
    return(assoc_out)
  }
  # Calculations only between nominal and numeric variables
  if (method == "eta_squared") {
    is_numeric <- which(sapply(data, class) %in% c("integer", "numeric"))
    if (length(is_numeric) == 0) {
      stop("Error: no numeric variables in `data`")
    }
    numeric_vars <- colnames(data)[is_numeric]
    is_factor <- which(sapply(data, class) == "factor")
    if (length(is_factor) == 0) {
      stop("Error: no factor variables in `data`")
    }
    factor_vars <- colnames(data)[is_factor]
    anova_pairs <- expand.grid(numeric_var = numeric_vars, factor_var = factor_vars)
    eta_squared_vec <- apply(
      anova_pairs, 
      1, 
      function(x) {
        formula_str <- paste0(x[1], "~", x[2])
        aov_fit <- aov(as.formula(formula_str), data = data)
        return(lsr::etaSquared(aov_fit)[, "eta.sq"])
      },
      simplify = T
    )
    eta_squared_df <- cbind(anova_pairs, eta_squared = eta_squared_vec)
    # Reshape into a matrix
    assoc_out <- matrix(NA, nrow = length(factor_vars), ncol = length(numeric_vars))
    rownames(assoc_out) <- factor_vars
    colnames(assoc_out) <- numeric_vars
    for (i in factor_vars) {
      for (j in numeric_vars) {
        row_index <- (eta_squared_df$factor_var == i & eta_squared_df$numeric_var == j)
        assoc_out[i, j] <- eta_squared_df[row_index, "eta_squared"]
      }
    }
    return(assoc_out)
  }
  #TODO: Cramer's V calculations
  # Calculations only between nominal variables
  if (method == "cramers_v") {
    is_factor <- sapply(data, class) == "factor"
  }
}
