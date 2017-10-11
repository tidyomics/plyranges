# io-bed.R

#' Read a BED file
#'
#' @param file A path to a file or a connection.
#' @param col_names An optional character vector for including additional
#' columns in \code{file} that are not part of the BED specification.
#' @param genome_info An optional character string or a Ranges object
#' that contains information about the genome build. For example the USSC identifier
#'"hg19" will add build information to the returned GRanges.
#' @param overlap_ranges An optional Ranges object. Only the intervals in the file
#' that overlap the Ranges will be returned.
#'
#' @description This is a lightweight wrapper to the import family
#' of functions defined in \pkg{rtracklayer}.
#'
#' @return A GRanges object
#'
#' @importFrom rtracklayer import.bed
#' @importFrom GenomeInfoDb seqinfo
#' @seealso \code{\link[rtracklayer]{BEDFile}}
#' @export
#' @rdname bed-files-read
read_bed <- function(file, col_names = NULL, genome_info = NULL,
                     overlap_ranges = NULL) {
  if (is.null(genome_info)) { genome_info <- NA }
  if (is(genome_info, "GRanges")) {
    seq_info <- seqinfo(genome_info)
    genome_info <- NA
  } else {
    seq_info <- NULL
  }
  import.bed(file, colnames = col_names,
             genome = genome_info,
             seqinfo = seq_info,
             which = overlap_ranges)
}


#' Write a BED file
#'
#' @param x A GRanges object
#' @param path Path or connection to write to
#'
#' @description This is a lightweight wrapper to the export family
#' of functions defined in \pkg{rtracklayer}.
#'
#' @importFrom rtracklayer export.bed
#' @seealso \code{\link[rtracklayer]{BEDFile}}
#' @export
#' @rdname bed-files-write
write_bed <- function(x, path) {
  export.bed(x, path)
}


