% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_utils.R
\name{ggpca}
\alias{ggpca}
\title{ggplot2-based PCA plot}
\usage{
ggpca(
  data,
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
  return_data = F
)
}
\arguments{
\item{data}{A SummarizedExperiment-like object. Must be compatible with \code{assay()} and
\code{colData()}. The columns in \code{assay} equate to data points in the PC scatterplot.}

\item{group_var}{A string or vector of strings for the grouping variable(s) to use from
\code{colData} to color points. If multiple variables are specified, they are combined into a
single string to make a single new variable.}

\item{pc_x}{A numeric. The PC to plot on the x-axis.}

\item{pc_y}{A numeric. The PC to plot on the y-axis.}

\item{ntop}{A numeric. Specifies the top \code{ntop} rows ranked by decreasing variance to subset the
data to prior to PCA.}

\item{center}{A logical. Should the data be zero-centered prior to PCA?}

\item{scale}{A logical. Should the data be scaled to unit variance prior to PCA?}

\item{equal_axes}{A logical. Should the scatterplot use the same axis limits for both axes?}

\item{point_size}{A numeric for the plot point size.}

\item{alpha}{A numeric for the plot point alpha level.}

\item{color}{A string for the color of plot points. Only used when \code{group_var} is \code{NULL}.}

\item{palette}{A string for the RColorBrewer palette name to use when \code{group_var} is specified.}

\item{return_data}{A logical. Should plot data be returned instead of a ggplot object?}
}
\value{
A ggplot object unless \code{return_data} is \code{TRUE}, then a data frame with the
user-specified PCs, grouping variable, and an attribute for the percent variace explained for
each user-specified PC.
}
\description{
Create an eigenvectors scatterplot from principal component analysis.
}
\details{
Conducts a principal component analysis (PCA) then produces a scatterplot using the
user-specified eigenvectors (ie., principal components or PCs) from the eigenvector matrix. Data
points can be optionally color-coded based on a user-specified variable.
}
\seealso{
\code{\link{prcomp}} \code{\link{SummarizedExperiment}} \code{\link{assay}}
\code{\link{colData}}
}
