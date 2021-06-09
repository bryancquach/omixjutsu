# Collection of functions for GENCODE data munging
# Author: Bryan Quach (bquach@rti.org)

#' Subset GENCODE GTF to specific features
#'
#' Creates a data frame with a subset of GTF records based on feature type.
#'
#' Filters a GTF-format data frame to only records of a specific feature type (based on column 3
#' of a GTF record). Adds an additional column with the GENCODE ID parsed from column 9.
#' 
#' @param gtf A data frame with data from an imported GENCODE GTF file.
#' @param feature_type A string denoting the feature type for which to retrieve records. This 
#'   corresponds to the 3rd column values in a GTF file.
#' @return A data frame in GTF-like format with an additional column (column 10) for the 
#'   GENCODE ID.
#' @export
subset_gencode_gtf <- function(gtf, feature_type = c("transcript", "gene", "exon")){
  feature_type <- match.arg(feature_type)
  gtf_subset <- gtf[gtf[, 3, drop = T] == feature_type, ]
  gtf_subset[, 10] <- switch (
    feature_type,
    "transcript" = {
      gsub(
        sapply(
          strsplit(gtf_subset[, 9, drop = T], split = ";", fixed = T), 
          function(split_data){
            split_data[[1]][1]
          },
          simplify = T
        ), 
        pattern = "transcript_id ", 
        replacement = ""
      )
    },
    "gene" = {
      gsub(
        sapply(
          strsplit(gtf_subset[, 9, drop = T], split = ";", fixed = T), 
          function(split_data){
            split_data[[1]][1]
          },
          simplify = T
        ), 
        pattern = "gene_id ", 
        replacement = ""
      )
    },
    "exon" = {
      gsub(
        sapply(
          strsplit(gtf_subset[, 9, drop = T], split = ";", fixed = T), 
          function(split_data){
            split_data[[1]][1]
          },
          simplify = T
        ), 
        pattern = "exon_id ", 
        replacement = ""
      )
    }
  )
  return(gtf_subset)
}

#' Retrieve mitochondrial feature IDs
#' 
#' Gets mitochondrial feature IDs from a GENCODE GTF-format data frame.
#' 
#' @inheritParams subset_gencode_gtf
#' @return A string vector with GENCODE IDs
#' @export
get_gencode_mitochondrial_ids <- function(gtf, feature_type = c("transcript", "gene", "exon")){
  gtf_subset <- subset_gencode_gtf(gtf, feature_type)
  ids <- gtf_subset[gtf_subset[, 1, drop = T] == "chrM", 10]
  return(ids)
}

#' Retrieve ribosomal RNA feature IDs
#' 
#' Gets ribosomal RNA feature IDs from a GENCODE GTF-format data frame.
#' 
#' @inheritParams subset_gencode_gtf
#' @return A string vector with GENCODE IDs
#' @export
get_gencode_rrna_ids <- function(gtf, feature_type = c("transcript", "gene", "exon")){
  gtf_subset <- subset_gencode_gtf(gtf, feature_type)
  ribo_index <- grep(
    sapply(
      strsplit(gtf_subset$V9, split = ";", fixed = T), 
      function(x){
        x[[2]][1]
      }, 
      simplify = T
    ), 
    pattern = " rRNA$", 
    perl = T
  )
  ids <- gtf_subset[ribo_index, 10]
  return(ids)
}

#' Retrieve protein-coding chrY gene IDs
#' 
#' Gets protein-coding chrY gene IDs from a GENCODE GTF-format data frame.
#' 
#' @param gtf A data frame with data from an imported GENCODE GTF file.
#' @param exclude_par A logical. Should genes in the pseudo-autosomal region be excluded?
#' @return A string vector with GENCODE IDs
#' @export
get_gencode_chr_y_gene_ids <- function(gtf, exclude_par = F){
  gtf_subset <- subset_gencode_gtf(gtf, "gene")
  chr_y_index <- which(gtf_subset[, 1, drop = T] == "chrY")
  protein_coding_index <- which(
    grepl(gtf_subset[, 9, drop = T], pattern = "gene_type protein_coding", fixed = T)
  )
  chr_y_protein_coding_index <- intersect(chr_y_index, protein_coding_index)
  chr_y_ids <- gtf_subset[chr_y_protein_coding_index, 10, drop = T]
  if (exclude_par) {
    chr_y_ids <- chr_y_ids[grep(chr_y_ids, pattern = "PAR_Y", fixed = T, invert = T)]
  }
  return(chr_y_ids)
}

#' Create GENCODE ID to gene name mapping.
#'
#' Creates a GENCODE ID to gene name mapping using a GTF file.
#'
#' @param gtf_file A string denoting a GTF file from which to build the mapping.
#' @param output_file Optional path and name of file to which to save the mapping table.
#' @return A data frame with each row linking a GENCODE ID to a gene name and HGNC ID.
#' @export
make_gencode_gene_id_map <- function(gtf_file, output_file = NULL){
  if (! file.exists(gtf_file)) {
    stop("Error: GTF file not found")
  }
  cat("\nLoading GTF file...\n")
  gtf <- read.delim(gtf_file, comment.char = "#'", sep = "\t", header = F, stringsAsFactors = F)
  cat("\nParsing GTF file...\n")
  gene_row_index <- (gtf[, 3] == "gene")
  gene_gtf_col9 <- strsplit(gtf[gene_row_index, 9], split = ";")
  mapping_list <- lapply(
    1:length(gene_gtf_col9),
    function(i){
      col9_tmp <- gene_gtf_col9[[i]]
      gene_id_data <- col9_tmp[grep(col9_tmp, pattern = "gene_id")]
      gene_id_data <- strsplit(gene_id_data, split = " ")[[1]]
      gene_id <- tail(gene_id_data, n = 1)
      gene_name_data <- col9_tmp[grep(col9_tmp, pattern = "gene_name")]
      gene_name_data <- strsplit(gene_name_data, split = " ")[[1]]
      gene_name <- tail(gene_name_data, n = 1)
      hgnc_id_data <- col9_tmp[grep(col9_tmp, pattern = "hgnc_id")]
      if (length(hgnc_id_data) == 0) {
        hgnc_id <- NA
      } else {
        hgnc_id_data <- strsplit(hgnc_id_data, split = " ")[[1]]
        hgnc_id <- tail(hgnc_id_data, n = 1)
      }
      return(c(gene_id, hgnc_id, gene_name))
    }
  )
  mapping_table <- as.data.frame(do.call(rbind, unique(mapping_list)), stringsAsFactors = F)
  colnames(mapping_table) <- c("gencode_id", "hgnc_id", "gene_name")
  if (! is.null(output_file)) {
    cat("\nWriting mapping table to file...\n")
    write.table(mapping_table, output_file, sep = "\t", quote = F, col.names = T, row.names = F)
  }
  return(mapping_table)
}
