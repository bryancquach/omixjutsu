% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnaseq_workflow_qc.R
\name{plot_dna_contamination_ratio}
\alias{plot_dna_contamination_ratio}
\title{Plot DNA contamination ratio}
\usage{
plot_dna_contamination_ratio(
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
  x_title = "Intergenic-intronic mappping ratio",
  y_title = "Sample count"
)
}
\arguments{
\item{data}{A two column data frame with intergenic mapping rates as the first column and
intronic mapping rates as the second column.}

\item{group}{A single column data frame denoting categorical variable values to use to
partition the data points into separate distributions/boxplots. Row names should match
rownames in \code{data}.}

\item{ids}{A vector of row names for subsetting \code{data} for plotting.}

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
A list with two ggplot objects \code{hist} and \code{boxplot} corresponding to a histogram and
jittered boxplot respectively. If \code{return_data} is \code{TRUE}, then a data frame with the plot
data values and grouping variable.
}
\description{
Create an intergenic-intronic mapping proportion histogram and boxplot.
}
\details{
Calculates the ratio between intergenic and intronic mapping rates and plots this as a histogram
and boxplot. The two plots are intended to be used as a single column, two row figure panel. The
intergenic and intronic mapping rates can be calculated uing
\code{\link{plot_mapping_categories}}.
}
\seealso{
\code{\link{plot_mapping_categories}} \code{\link{hist_boxplot2}}
}
