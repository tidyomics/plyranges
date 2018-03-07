#' Read a BigWig file
#' @param file A path to a file or URL.
#' @param genome_info An optional character string or a Ranges object
#' that contains information about the genome build. For example the identifier
#'"hg19" will add build information to the returned GRanges.
#' @param overlap_ranges An optional Ranges object. Only the intervals in the file
#' that overlap the Ranges will be loaded.
#'
#' @examples
#' if (.Platform$OS.type != "windows") {
#'   test_path <- system.file("tests", package = "rtracklayer")
#'   bw_file <- file.path(test_path, "test.bw")
#'   gr <- read_bigwig(bw_file)
#'   gr
#' }
#' @return a GRanges object
#' @importFrom rtracklayer import.bw BigWigSelection BigWigFile
#' @seealso \link[rtracklayer]{BigWigFile}
#' @export
#' @rdname io-bigwig-read
read_bigwig <- function(file, genome_info = NULL, overlap_ranges = NULL) {

  stopifnot(is(overlap_ranges, "GRanges") || is.null(overlap_ranges))

  if (!is.null(genome_info)) {
    if (is(genome_info, "GRanges")) {
      seq_info <- seqinfo(genome_info)
    } else if (is.character(genome_info)) {
      seq_info <- GenomeInfoDb::Seqinfo(genome = genome_info)
    }
  } else {
      seq_info <- NULL
  }

  if (!is.null(overlap_ranges)) {
    selection <- BigWigSelection(overlap_ranges)
    rng <- import.bw(file, selection = selection, as = "GRanges")
  } else {
    rng <- import.bw(file, format = "bw", as = "GRanges")
  }

  if (!is.null(seq_info)) {
    seqinfo(rng) <- seq_info
  }
  rng
}

#' Write a BED file
#'
#' @param x A GRanges object
#' @param file File name, URL or connection specifying a file to write x to.
#'             Compressed files with extensions such as '.gz' are handled
#'             automatically.
#'
#' @description This is a lightweight wrapper to the export family
#' of functions defined in \pkg{rtracklayer}.
#'
#' @importFrom rtracklayer export.bw
#' @seealso \code{\link[rtracklayer]{BigWigFile}}
#' @export
#' @return The write functions return a BigWigFile invisibly
#' @rdname io-bigwig-write
write_bigwig <- function(x, file) {
  export.bw(x, file)
}
