% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnaseq_workflow_qc.R
\name{plot_hisat2_stats}
\alias{plot_hisat2_stats}
\title{Plot HISAT2 mapping statistics as a scatterplot}
\usage{
plot_hisat2_stats(
  data,
  point_size = 4,
  fill = "red4",
  alpha = 0.75,
  ids = NULL,
  invert = F
)
}
\arguments{
\item{data}{A data frame HISAT2 mapping data in the same format as \code{hisat2} returned by
\code{\link{load_paired_end_qc_data}}.}

\item{point_size}{A numeric. The size for data points.}

\item{fill}{A string. The fill color for data points.}

\item{alpha}{A numeric. The alpha level for data points.}

\item{ids}{A vector of rownames for subsetting \code{data} for plotting.}

\item{invert}{A logical. Should \code{ids} be used for excluding rows from plotting instead?}
}
\value{
A ggplot object.
}
\description{
Create a concordant alignment vs. overall alignment rate scatterplot. Only applicable for
paired-end sequencing data.
}
\seealso{
\code{\link{load_paired_end_qc_data}}
}
