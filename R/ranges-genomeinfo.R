# Representing seqinfo information as a Ranges

#' Construct annotation information for Ranges
#'
#' @param genome A character vector of length one indicating the genome build.
#' If this is the only argument supplied, the build information will be
#' retrieved from UCSC database.
#' @param seqnames A character vector containing the name of sequences.
#' @param seqlengths An optional integer vector containg the lengths of sequences.
#' @param is_circular An optional logical vector indicating whether a sequence is ciruclar.
#'
#' @return a GRanges object
#' @importFrom GenomeInfoDb Seqinfo
#' @export
genome_info <- function(genome = NULL, seqnames = NULL, seqlengths = NULL, is_circular = NULL) {

  seqlengths <- ifelse(length(seqlengths) == 0, NA, seqlengths)
  is_circular <- ifelse(length(is_circular) == 0, NA, is_circular)
  genome <- ifelse(length(genome) == 0, NA, genome)
  seqinfo_res <- Seqinfo(seqnames, seqlengths, is_circular, genome)
  get_genome_info(seqinfo_res)

}

#' Add annotation information to a Ranges object
#'
#' @param .data
#' @param genome A character vector of length one indicating the genome build.
#' @param seqnames A character vector containing the name of sequences.
#' @param seqlengths An optional integer vector containg the lengths of sequences.
#' @param is_circular An optional logical vector indicating whether a sequence is ciruclar.
#'
#' @return An annotated GRanges object. To see the annotations use
#' \code{get_genome_info}
#' @export
set_genome_info <- function(.data, genome = NULL, seqnames = NULL,
                            seqlengths = NULL, is_circular = NULL) {

  info <- seqinfo(genome_info(genome, seqnames, seqlengths, is_circular))
  seqinfo(.data) <- info
}


#' Extract annotation as a Ranges object
#'
#' @param .data a Seqinfo object or a GRanges with a seqinfo slot
#'
#' @return A GRanges object
#' @export
get_genome_info <- function(.data)  UseMethod("get_genome_info")

get_genome_info.GenomicRanges <- function(.data) {
  info <- seqinfo(.data)
  get_genome_info(info)
}

get_genome_info.Seqinfo <- function(.data) {

  if (any(is.na(seqlengths(.data)))) {
    stop("seqlengths must be non-missing to convert to GRanges")
  }
  as(.data, "GRanges")
}

