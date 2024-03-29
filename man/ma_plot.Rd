% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deseq2.R, R/dge_utils.R
\name{ma_plot.DESeqResults}
\alias{ma_plot.DESeqResults}
\alias{ma_plot}
\alias{ma_plot.data.frame}
\title{Bland-Altman plot}
\usage{
\method{ma_plot}{DESeqResults}(
  dds_results,
  outliers = NULL,
  sig_cutoff = 0.05,
  sig_up_color = "red4",
  sig_down_color = "steelblue3",
  nonsig_color = "gray60",
  outlier_color = "goldenrod3",
  sig_alpha = 0.8,
  nonsig_alpha = 0.4,
  sig_size = 3,
  nonsig_size = 2,
  outlier_size = sig_size * 2.5
)

ma_plot(object, ...)

\method{ma_plot}{data.frame}(
  data,
  basemean_colname,
  fc_colname,
  adj_pvalue_colname,
  outlier_index = NULL,
  sig_cutoff = 0.05,
  sig_up_color = "red4",
  sig_down_color = "steelblue3",
  nonsig_color = "gray60",
  outlier_color = "goldenrod3",
  sig_alpha = 0.8,
  nonsig_alpha = 0.4,
  sig_size = 3,
  nonsig_size = 2,
  outlier_size = sig_size * 2.5
)
}
\arguments{
\item{outliers}{A vector of outlier IDs.}

\item{sig_cutoff}{Adjusted p-value significance threshold.}

\item{sig_up_color}{Point color for statistically significant up-regulated features.}

\item{sig_down_color}{Point color for statistically significant down-regulated features.}

\item{nonsig_color}{Point color for features not statistical significant.}

\item{outlier_color}{Point color for statistically significant features flagged as outliers.}

\item{sig_alpha}{Point color alpha value for statistically significant features.}

\item{nonsig_alpha}{Point color alpha value for features not statistically significant.}

\item{sig_size}{Point size value for statistically significant features.}

\item{nonsig_size}{Point size value for outliers.}

\item{object}{A DESeqResults or data.frame object.}

\item{basemean_colname}{A string denoting the column in \code{data} that contains the base mean
expression values.}

\item{fc_colname}{A string denoting the column in \code{data} that contains the log2 fold changes.}

\item{adj_pvalue_colname}{A string denoting the column in \code{data} that contains the adjusted
p-values.}

\item{outlier_index}{A numeric vector denoting row indices of outliers.}
}
\value{
A ggplot2 object.
}
\description{
Plot log fold changes vs. log expression level.
}
\details{
Plots the log2 fold change vs log10 expression level for each feature tested for
differential expression.
}
