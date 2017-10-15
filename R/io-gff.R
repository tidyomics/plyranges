#' Read a GFF/GTF/GVT file
#'
#' @param file A path to a file or a connection.
#' @param col_names An optional character vector for parsing specific
#' columns in \code{file} that are part of the GFF specification.
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
#' @importFrom rtracklayer import.gff
#' @importFrom GenomeInfoDb seqinfo
#' @seealso \code{\link[rtracklayer]{GFFFile}}
#' @export
#' @rdname gff-files-read
#' @importFrom rtracklayer import.gff
read_gff <- function(file, col_names = NULL, genome_info = NULL, overlap_ranges = NULL) {
  if (is(genome_info, "GRanges")) {
    genome_info <- seqinfo(genome_info)
  } else {
    genome_info <- NA
  }

  import.gff(file,
             colnames = col_names,
             genome = genome_info,
             which = overlap_ranges)
}

#' @export
#' @rdname gff-files-read
#' @importFrom rtracklayer import.gff1
read_gff1 <- function(file, col_names = NULL, genome_info = NULL, overlap_ranges = NULL) {
  if (is(genome_info, "GRanges")) {
    genome_info <- seqinfo(genome_info)
  } else {
    genome_info <- NA
  }

  import.gff1(file,
              colnames = col_names,
              genome = genome_info,
              which = overlap_ranges)
}

#' @export
#' @rdname gff-files-read
#' @importFrom rtracklayer import.gff2
read_gff2 <- function(file, col_names = NULL, genome_info = NULL, overlap_ranges = NULL) {
  if (is(genome_info, "GRanges")) {
    genome_info <- seqinfo(genome_info)
  } else {
    genome_info <- NA
  }

  import.gff2(file,
              colnames = col_names,
              genome = genome_info,
              which = overlap_ranges)
}

#' @export
#' @rdname gff-files-read
#' @importFrom rtracklayer import.gff3
read_gff3 <- function(file, col_names = NULL, genome_info = NULL, overlap_ranges = NULL) {
  if (is(genome_info, "GRanges")) {
    genome_info <- seqinfo(genome_info)
  } else {
    genome_info <- NA
  }

  import.gff3(file,
              colnames = col_names,
              genome = genome_info,
              which = overlap_ranges)
}



#' Write a GFF(123) file
#'
#' @param x A GRanges object
#' @param path Path or connection to write to
#' @param index If TRUE the output file will be compressed and indexed using bgzf and
#' tabix.
#'
#' @description This is a lightweight wrapper to the export family
#' of functions defined in \pkg{rtracklayer}.
#'
#' @importFrom rtracklayer export.gff export.gff1 export.gff2 export.gff3
#' @seealso \code{\link[rtracklayer]{GFFFile}}
#' @export
#' @rdname gff-files-write
write_gff <- function(x, path, index = FALSE) {
  export.gff(x, path, index = index)
}

#' @export
#' @rdname gff-files-write
write_gff1 <- function(x, path, index = FALSE) {
  export.gff1(x, path, index = index)
}

write_gff2 <- function(x, path, index = FALSE) {
  export.gff2(x, path, index = index)
}

write_gff3 <- function(x, path, index = FALSE) {
  export.gff3(x, path, index = index)
}
