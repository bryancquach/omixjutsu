% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/table_utils.R
\name{select_replicates}
\alias{select_replicates}
\title{Choose best replicates to keep}
\usage{
select_replicates(data, id, metric, max_replicates = 1, decreasing = T)
}
\arguments{
\item{data}{Data frame with at least two columns, one with a sample ID for detecting replicates
and another with a metric to use for selecting replicates.}

\item{id}{Sample ID to use for detecting replicates.}

\item{metric}{Metric name to use for selecting replicates.}

\item{max_replicates}{Maximum number of replicates to keep from a replicate set.}

\item{decreasing}{If 'TRUE', choose replicates by starting with the highest metric value.
Otherwise choose replicates by starting with the lowest metric value.}
}
\value{
A vector of indices for selected sample replicates as well as samples without
replicates. Indices correspond to rows in \code{data}.
}
\description{
Reduce sample list by subsetting replicates based on a specified metric.
}
\details{
Detects replicate samples from a sample list using ID duplicates and retains only replicate(s)
that meet selection criteria for the user-specified metric.
}
