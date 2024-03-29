% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/table_utils.R
\name{below_max_cutoff}
\alias{below_max_cutoff}
\title{Detect samples below a specified cutoff value}
\usage{
below_max_cutoff(
  metrics_table,
  metric,
  cutoff = NULL,
  include_equals = F,
  na_as_false = F,
  coef = 1.7,
  use_unskewed_iqr = F,
  use_skewed_iqr = F
)
}
\arguments{
\item{metrics_table}{Data frame or matrix with metrics as columns and samples as rows.}

\item{metric}{Name of column with metric to use.}

\item{cutoff}{Cutoff value.}

\item{include_equals}{If \code{TRUE}, designate metric values equal to cutoff value as TRUE.}

\item{na_as_false}{If \code{TRUE}, designate comparisons to 'NA' values as FALSE.}

\item{coef}{Constant to multiply the interquartile range by in calculating a relative cutoff.
The default value of 1.7 approximately equates to values at least 3 standard deviations from
the mean for a gaussian distribution.}

\item{use_unskewed_iqr}{If \code{TRUE}, set a relative cutoff value using the interquartile range
of the metric distribution multiplied by \code{coef}. More suitable for non-skewed distributions.}

\item{use_skewed_iqr}{If \code{TRUE}, set relative cutoff value using a medcouple-based adjustment
of the interquartile range of the metric distribution multipled by \code{coef}. More suitable for
skewed distributions.}
}
\value{
A list containing the following elements:
\itemize{
\item \code{keep} a logical vector indicating whether a given value for a metric meets the designated
cutoff criteria. The vector order matches the input row order of \code{metrics_table}.
\item \code{cutoff} The cutoff value used.
}
}
\description{
Identify samples with values below a maximum cutoff value.
}
\details{
Indicates whether a sample is below a maximum cutoff value for a user-specified metric. Includes
options to set the cutoff using values derived from interquartile range for both skewed and
unskewed distributions.
}
