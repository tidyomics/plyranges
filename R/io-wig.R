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
#' @seealso [rtracklayer::WIGFile()]
#' @rdname io-wig-read
#' @examples
#' test_path <- system.file("tests", package = "rtracklayer")
#' test_wig <- file.path(test_path, "step.wig")
#' gr <- read_wig(test_wig)
#' gr
#' gr <- read_wig(test_wig, genome_info = "hg19")
#' @return A GRanges object
#' @export
read_wig <- function(file, genome_info = NULL, overlap_ranges = NULL) {

  args <- norm_args_reader(genome_info)

  import.wig(file,
            trackLine = FALSE,
            genome = args$genome_info,
            seqinfo = args$seq_info,
            which = overlap_ranges)
}


#' Write a WIG file
#'
#' @param x A GRanges object
#' @param file File name, URL or connection specifying a file to write x to.
#'             Compressed files with extensions such as '.gz' are handled
#'             automatically.
#' @importFrom rtracklayer export.wig
#' @seealso [rtracklayer::WIGFile()]
#' @export
#' @return The write function returns a WIGFile invisibly.
#' @rdname io-wig-write
write_wig <- function(x, file) {
  export.wig(x, file)
}
