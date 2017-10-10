#' Read a BEDGraph file
#' @rdname bed-files-read
#' @export
#' @importFrom rtracklayer import.bedGraph
read_bed_graph <- function(file, col_names = NULL, genome_info = NULL,
                           overlap_ranges = NULL) {
  if (is.null(genome_info)) { genome_info <- NA }
  import.bedGraph(file, colnames = col_names,
             genome = genome_info,
             which = overlap_ranges)
}

#' Write a BEDGraph file
#'
#' @importFrom rtracklayer export.bedGraph
#' @export
#' @rdname bed-files-write
write_bed <- function(x, path) {
  export.bed(x, path)
}
