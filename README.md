# R package: `omixjutsu`

The `omixjutsu` R package is a suite of utility functions for munging, analyzing, and visualizing various omics data. The following types of tasks are currently supported (more to come!):

* Basic data visualization (histograms, boxplots, barplots, scatterplots, heatmpas, etc.) using a ggplot2 framework.
* Munging and visualization of RNA-seq data processing (read trimming, alignment, transcript quantification, etc.) QC metrics as generated by [MultiQC](https://multiqc.info/).
* Differential gene expression analysis using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).
* Differential gene expression analysis using [limma and voom](https://bioconductor.org/packages/limma/).
* Evaluating pairwise variable relationships using correlation metrics, Eta-squared, and Cramer's V.
* Evaluating linear, additive model fit using percentage of explained gene expresion variance (PVE analysis).
* [Independent hypothesis weighting](https://bioconductor.org/packages/IHW/) for multiple testing correction.

# Installation

The most recent development version of this package can be installed using [`devtools`](https://devtools.r-lib.org/) as follows:

```R
devtools::install_github("bryancquach/omixjutsu")
```

To install a version from a specific commit or branch use the `ref` option:

```R
# Replace <id> with the desired commit ID (i.e., the SHA hash) or repository branch
devtools::install_github("bryancquach/omixjutsu", ref = "<id>")

# Example for a branch named 'example_branch'
devtools::install_github("bryancquach/omixjutsu", ref = "example_branch")
```

See [this GitHub issue](https://github.com/bryancquach/omixjutsu/issues/2#issue-1557056581) for troubleshooting common sources of installation failures.
