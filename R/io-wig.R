#' Read a WIG file
#'
#' @param file A path to a file or a connection.
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
#' @importFrom rtracklayer import.wig
#' @seealso \code{\link[rtracklayer]{BEDFile}}
#' @export
#' @rdname bed-files-read
#' @importFrom rtracklayer import.wig
#' @importFrom GenomeInfoDb seqinfo
read_wig <- function(file,  genome_info = NULL, overlap_ranges = NULL) {

  if (is.null(genome_info)) { genome_info <- NA }
  if (is(genome_info, "GRanges")) {
    seq_info <- seqinfo(genome_info)
    genome_info <- NA
  } else {
    seq_info <- NULL
  }

  import.wig(file,
            trackLine = FALSE,
            genome = genome_info,
            seqinfo = seq_info,
            which = overlap_ranges)
}
