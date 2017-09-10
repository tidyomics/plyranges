# readr read_* equivalents
# rtracklayer and Rsamtools cleaner API, they always return
# a GRanges object, perhaps which should be changed to overlap?
# read_bam

# read_bw
#' @importFrom rtracklayer import.bw
read_bigwig <- function(file, ...) {
  import.bw(file, ...)
}

#' @importFrom rtracklayer import.wig
read_wig <- function(file,  genome = NULL, seqinfo = NULL, which = NULL) {
  import.wig(file, genome = genome, seqinfo = NULL, which = NULL)
}

#' Read a BED file
#'
#' @param file A path to a file or a connection.
#' @param col_names An optional character vector for including additional
#' columns in \code{file} that are not part of the BED specification.
#' @param genome_info An optional character string or a Ranges object
#' that contains information about the genome build. For example the USSC identifier
#'"hg19".
#' @param overlap_ranges An optional Ranges object. Only the intervals in the file
#' that overlap the Ranges will be returned.
#'
#' @description This is a lightweight wrapper to the import family
#' of functions defined in \pkg{rtracklayer}.
#'
#' @importFrom rtracklayer import.bed
#' @seealso \code{\link[rtracklayer]{BEDFile}}
#' @export
read_bed <- function(file, col_names = NULL, genome_info = NULL,
                     overlap_ranges = NULL) {
  if (is.null(genome_info)) { genome_info <- NA }
  import.bed(file, colnames = col_names,
             genome = genome_info,
             which = overlap_ranges)
}

#' @importFrom rtracklayer import.bedGraph
read_bed_graph <- function(file, col_names = NULL, genome = NULL, seqinfo = NULL,
                           which = NULL, col_extra = NULL) {
  import.bedGraph(file, genome = genome, colnames = col_names, genome = genome,
             seqinfo = seqinfo, which = which, extraCols = col_extra)
}

#' @importFrom rtracklayer import.gff
read_gff <- function(file, col_names = NULL, genome = NULL, which = NULL, feature_type = NULL) {
  import.gff(file, colnames = col_names, genome = genome, which = which,
             feature_type = feature_type)
}

#' @importFrom rtracklayer import.gff1
read_gff1 <- function(file, col_names = NULL, genome = NULL, which = NULL, feature_type = NULL) {
  import.gff1(file, colnames = col_names, genome = genome, which = which,
              feature_type = feature_type)
}

#' @importFrom rtracklayer import.gff2
read_gff2 <- function(file, col_names = NULL, genome = NULL, which = NULL, feature_type = NULL) {
  import.gff2(file, colnames = col_names, genome = genome, which = which,
              feature_type = feature_type)
}

#' @importFrom rtracklayer import.gff3
read_gff3 <- function(file, col_names = NULL, genome = NULL, which = NULL, feature_type = NULL) {
  import.gff3(file, colnames = col_names, genome = genome, which = which,
              feature_type = feature_type)
}

