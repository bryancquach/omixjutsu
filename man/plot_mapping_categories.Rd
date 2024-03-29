% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnaseq_workflow_qc.R
\name{plot_mapping_categories}
\alias{plot_mapping_categories}
\title{Plot annotation category percentages for alignments}
\usage{
plot_mapping_categories(
  data,
  fill = NULL,
  sort = F,
  ids = NULL,
  invert = F,
  consolidate = F,
  return_data = F
)
}
\arguments{
\item{data}{A data frame with annotation data in the same format as \code{rseqc_alignment_category}
returned by \code{\link{load_paired_end_qc_data}}.}

\item{fill}{A string vector with 3 or 6 colors corresponding to cateogires for the barplot.}

\item{sort}{A logical. Should the x-axis be ordered by coding sequence/exonic alignment
percentage?}

\item{ids}{A vector of rownames for subsetting \code{data} for plotting.}

\item{invert}{A logical. Should \code{ids} be used for excluding rows from plotting instead?}

\item{consolidate}{A logical. Should mapping categories be consolidated into the three more
general categories of exonic, intronic, and intergenic?}

\item{return_data}{A logical. Should plot data be returned instead of a ggplot object?}
}
\value{
A ggplot object.  If \code{return_data} is \code{TRUE}, then a data frame with the plot data
values for each genomic annotation category as columns.
}
\description{
Create stacked barplots with percentages for types of genomic annotations.
}
\details{
Create barplots with the percentage of alignments that map to a genomic region with annotations
of intergenic, proximally down and upstream of genes, intronic, UTRs or coding sequence.
Includes the option to consolidate categories into more general categories for intergenic,
intronic, and exonic regions.
}
\seealso{
\code{\link{load_paired_end_qc_data}}
}
