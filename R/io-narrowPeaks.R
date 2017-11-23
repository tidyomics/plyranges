# narrowpeaks

#' Read a narrowPeaks file
#'
#' @param file A path to a file or a connection..
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
#' @seealso \code{\link[rtracklayer]{BEDFile}}
#' @export
#' @rdname bedxy-files-read
read_narrowpeaks <- function(file, genome_info = NULL,
                     overlap_ranges = NULL) {
  if (is.null(genome_info)) { genome_info <- NA }

  if (is(genome_info, "GRanges")) {
    seq_info <- seqinfo(genome_info)
    genome_info <- NA
  } else {
    seq_info <- NULL
  }

  extra_cols <- c(rep("numeric", 3), "integer")
  col_names <- c("signalValue", "pvalue", "qvalue", "peak")
  import.bed(file, colnames = col_names,
             extraCols = extra_cols,
             genome = genome_info,
             which = overlap_ranges)
}

#' Write a narrowPeaks file
#'
#' @param x A GRanges object
#' @param path Path or connection to write to
#' @param col_names A character vector indicating which columns
#' in the Ranges object represent columns in the narrowPeaks spec
#'
#' @description This is a lightweight wrapper to the import family
#' of functions defined in \pkg{rtracklayer}.
#'
#' @importFrom rtracklayer export.bed
#' @seealso \code{\link[rtracklayer]{BEDFile}}
#' @export
#' @rdname bedxy-files-write
write_narrowpeaks <- function(x, path, col_names) {
  export.bed(x, path, expNames = col_names)
}
