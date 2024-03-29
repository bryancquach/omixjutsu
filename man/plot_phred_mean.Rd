% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnaseq_workflow_qc.R
\name{plot_phred_mean}
\alias{plot_phred_mean}
\title{Plot per sample mean of per read mean PHRED scores}
\usage{
plot_phred_mean(
  data,
  group = NULL,
  ids = NULL,
  invert = F,
  return_data = F,
  binsize = NULL,
  colors = NULL,
  hist_alpha = 0.75,
  box_fill = "goldenrod",
  box_alpha = 0.5,
  box_lwd = 1,
  jitter_alpha = 0.75,
  jitter_size = 1.75,
  x_title = "Mean of mean PHRED score",
  y_title = "Sample count"
)
}
\arguments{
\item{data}{A data frame with data in the same format as \code{phred_seq} returned by
\code{\link{load_paired_end_qc_data}}.}

\item{group}{A single column data frame denoting categorical variable values to use to
partition the data points into separate distributions/boxplots. Row names should match rownames
in \code{data}.}

\item{ids}{A vector of rownames for subsetting \code{data} for plotting.}

\item{invert}{A logical. Should \code{ids} be used for excluding rows from plotting instead?}

\item{return_data}{A logical. Should plot data be returned instead of a ggplot object?}

\item{binsize}{A numeric value for the histogram bin widths.}

\item{colors}{A string vector. The histogram bar colors and jitter colors for each categorical
variable value.}

\item{hist_alpha}{A numeric value for the alpha level for histogram bars.}

\item{box_fill}{A string. The fill color for the boxplot.}

\item{box_alpha}{A numeric value for the alpha level for the boxplot.}

\item{box_lwd}{A numeric value for the boxplot line width.}

\item{jitter_alpha}{A numeric for the alpha level for jittered data points.}

\item{jitter_size}{A numeric for the point size for jittered data points.}

\item{x_title}{A string denoting the x-axis title. Only added to the boxplot.}

\item{y_title}{A string denoting the y-axis title. Only added to the histogram.}
}
\value{
A ggplot object unless \code{return_data} is \code{TRUE}, then a data frame with the
mean of mean PHRED scores and grouping variable.
}
\description{
Create a histogram and boxplot of per sample mean of per read mean PHRED scores. The two plots
are intended to be used as a single column, two row figure panel.
}
\seealso{
\code{\link{load_paired_end_qc_data}} \code{\link{hist_boxplot2}}
}
