#' Read a BigWig file
#' @param file A path to a file or a connection.
#' @param genome_info An optional character string or a Ranges object
#' that contains information about the genome build. For example the USSC identifier
#'"hg19" will add build information to the returned GRanges.
#' @param overlap_ranges An optional Ranges object. Only the intervals in the file
#' that overlap the Ranges will be loaded.
#' @importFrom rtracklayer import.bw BigWigSelection BigWigFile
#' @seealso \code{\link[rtracklayer]{BigWigFile}}
#' @export
#' @rdname bigwig-files-read
read_bigwig <- function(file, genome_info = NULL, overlap_ranges = NULL) {

  file <- BigWigFile(file)

  if (!is.null(overlap_ranges) & is(overlap_ranges, "GRanges")) {
    selection <- BigWigSelection(overlap_ranges)
  } else {
    selection <- BigWigSelection(file)
  }

  if (is(genome_info, "GRanges")) {
    seq_info <- seqinfo(genome_info)
  } else if (!is.null(genome_info)){
    seq_info <- Seqinfo(genome = genome_info)
  } else {
    seq_info <- NULL
  }

  rng <- import.bw(file, selection = selection)
  if (!is.null(seq_info)) {
    seqinfo(rng) <- seq_info
    return(rng)
  }

  rng
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
#' @seealso \code{\link[rtracklayer]{BigWigFile}}
#' @export
#' @rdname bigwig-files-write
write_bed <- function(x, path) {
  export.bed(x, path)
}
