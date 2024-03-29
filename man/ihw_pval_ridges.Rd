% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ihw.R
\name{ihw_pval_ridges}
\alias{ihw_pval_ridges}
\title{Stratified p-value density distributions}
\usage{
ihw_pval_ridges(
  ihw_data,
  num_quantiles = 10,
  fill_color = "gray70",
  fill_alpha = 1
)
}
\arguments{
\item{ihw_data}{An ihwResult object.}

\item{num_quantiles}{Number of strata into which the independent variable is divided.}

\item{fill_color}{Fill color for area under density distributions.}

\item{fill_alpha}{Fill color alpha value for area under density distributions.}
}
\value{
A ggplot object.
}
\description{
Plot p-value density distribution for IHW variable strata.
}
\details{
This plot partitions the p-values by quantiles of the IHW variable to generate multiple p-value
density distributions. For an informative IHW variable, the larger the quantile ID (strata) the
more skew there should be towards enrichment of smaller p-values.
}
\seealso{
Other ihw: 
\code{\link{ihw_diagnostic_plots}()},
\code{\link{ihw_partitioned_pval_weights}()},
\code{\link{ihw_pval_cdf}()},
\code{\link{ihw_pval_comparison}()},
\code{\link{ihw_pval_decision_boundary}()},
\code{\link{ihw_pval_scatter}()},
\code{\link{ihw_stratified_pval_hist}()}
}
\concept{ihw}
