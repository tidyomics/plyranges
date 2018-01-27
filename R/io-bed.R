# io-bed.R
norm_args_reader <- function(genome_info) {

  if (is(genome_info, "GenomicRanges")) {
    seq_info <- seqinfo(genome_info)
    genome_info <- NA
  } else if (is.character(genome_info)) {
    seq_info <- NULL
  } else {
    genome_info <- NA
    seq_info <- NULL
  }

  return(list(genome_info = genome_info, seq_info = seq_info))
}

#' Read a BED or BEDGraph file
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
#' @importFrom methods is
#' @seealso \code{\link[rtracklayer]{BEDFile}}
#' @export
#' @rdname io-bed-read
read_bed <- function(file, col_names = NULL, genome_info = NULL,
                     overlap_ranges = NULL) {
  # check genome_info input
  args <- norm_args_reader(genome_info)

  import.bed(file, colnames = col_names,
             genome = args$genome_info,
             seqinfo = args$seq_info,
             which = overlap_ranges,
             trackLine = FALSE)
}


#' Write a BED file
#'
#' @param x A GRanges object
#' @param path Path or connection to write to
#' @param index Compress and index the output file
#'              with bgzf and tabix (default = FALSE). Note that tabix indexing will sort the
#'              data by chromosome and start.
#'
#' @description This is a lightweight wrapper to the export family
#' of functions defined in \pkg{rtracklayer}.
#'
#' @importFrom rtracklayer export.bed
#' @seealso \code{\link[rtracklayer]{BEDFile}}
#' @export
#' @rdname io-bed-write
write_bed <- function(x, path, index = FALSE) {
  export.bed(x, path, index = index)
}

#' @rdname io-bed-read
#' @export
#' @importFrom rtracklayer import.bedGraph
read_bed_graph <- function(file, col_names = NULL, genome_info = NULL,
                           overlap_ranges = NULL) {

  args <- norm_args_reader(genome_info)

  import.bedGraph(file, colnames = col_names,
                  genome = args$genome_info,
                  seqinfo = args$seq_info,
                  which = overlap_ranges,
                  trackLine = FALSE)
}


#' @importFrom rtracklayer export.bedGraph
#' @export
#' @rdname io-bed-write
write_bed_graph <- function(x, path, index = FALSE) {
  export.bedGraph(x, path, index = index)
}



