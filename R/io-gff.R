
norm_args_reader_gff <- function(genome_info) {
  if (!is.null(genome_info)) {

    if (is(genome_info, "GRanges")) {
      return(seqinfo(genome_info))
    }
    return(genome_info)
  }
  NA
}
#' Read a GFF/GTF/GVT file
#'
#' @param file A path to a file or a connection.
#' @param col_names An optional character vector for parsing specific
#' columns in \code{file} that are part of the GFF specification. These should
#' name either fixed fields, like source or type, or, for GFF2 and GFF3,
#' any attribute.
#' @param genome_info An optional character string or a Ranges object
#' that contains information about the genome build. For example the UCSC identifier
#'"hg19" will add build information to the returned GRanges.
#' @param overlap_ranges An optional Ranges object. Only the intervals in the file
#' that overlap the Ranges will be returned.
#'
#' @description This is a lightweight wrapper to the import family
#' of functions defined in \pkg{rtracklayer}.
#'
#' @return A GRanges object
#'
#' @importFrom rtracklayer import.gff
#' @importFrom GenomeInfoDb seqinfo
#' @seealso \code{\link[rtracklayer]{GFFFile}}
#' @export
#' @rdname io-gff-read
#' @importFrom rtracklayer import.gff
read_gff <- function(file, col_names = NULL, genome_info = NULL, overlap_ranges = NULL) {
  import.gff(file,
             colnames = col_names,
             genome = norm_args_reader_gff(genome_info),
             which = overlap_ranges)
}

#' @export
#' @rdname io-gff-read
#' @importFrom rtracklayer import.gff1
read_gff1 <- function(file, col_names = NULL, genome_info = NULL, overlap_ranges = NULL) {
  import.gff1(file,
              colnames = col_names,
              genome = norm_args_reader_gff(genome_info),
              which = overlap_ranges)
}

#' @export
#' @rdname io-gff-read
#' @importFrom rtracklayer import.gff2
read_gff2 <- function(file, col_names = NULL, genome_info = NULL, overlap_ranges = NULL) {
  import.gff2(file,
              colnames = col_names,
              genome = norm_args_reader_gff(genome_info),
              which = overlap_ranges)
}

#' @export
#' @rdname io-gff-read
#' @importFrom rtracklayer import.gff3
read_gff3 <- function(file, col_names = NULL, genome_info = NULL, overlap_ranges = NULL) {
  import.gff3(file,
              colnames = col_names,
              genome = norm_args_reader_gff(genome_info),
              which = overlap_ranges)
}



#' Write a GFF(123) file
#'
#' @param x A GRanges object
#' @param file Path or connection to write to
#' @param index If TRUE the output file will be compressed and indexed using bgzf and
#' tabix.
#'
#' @description This is a lightweight wrapper to the export family
#' of functions defined in \pkg{rtracklayer}.
#'
#' @importFrom rtracklayer export.gff export.gff1 export.gff2 export.gff3
#' @seealso \code{\link[rtracklayer]{GFFFile}}
#' @export
#' @rdname io-gff-write
write_gff <- function(x, file, index = FALSE) {
  export.gff(x, file, index = index)
}

#' @export
#' @rdname io-gff-write
write_gff1 <- function(x, file, index = FALSE) {
  export.gff1(x, file, index = index)
}

#' @export
#' @rdname io-gff-write
write_gff2 <- function(x, file, index = FALSE) {
  export.gff2(x, file, index = index)
}

#' @export
#' @rdname io-gff-write
write_gff3 <- function(x, file, index = FALSE) {
  export.gff3(x, file, index = index)
}
