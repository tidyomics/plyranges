
#' Filter by overlapping/non-overlapping ranges
#'
#' @param x,y Objects representing ranges
#' @param maxgap,minoverlap The maximimum gap between intervals as an integer
#' greater than or equal to zero. The minimum amount of overlap between intervals
#' as an integer greater than zero, accounting for the maximum gap.
#'
#' @importFrom IRanges subsetByOverlaps
#' @export
#' @rdname filter-overlaps.rd
filter_by_overlaps <- function(x,y, maxgap, minoverlap) {
  UseMethod("filter_by_overlaps")
}

filter_by_overlaps.Ranges <- function(x,y, maxgap = 0L, minoverlap =1L) {
  subsetByOverlaps(x,y,maxgap, minoverlap, type = "any")
}

filter_by_overlaps.GenomicRanges <- function(x,y, maxgap = 0L, minoverlap =1L) {
  subsetByOverlaps(x,y,maxgap, minoverlap, type = "any", ignore.strand = TRUE)
}

#' @export
#' @rdname filter-overlaps.rd
filter_by_non_overlaps <- function(x,y, maxgap, minoverlap) {
  UseMethod("filter_by_non_overlaps")
}

filter_by_non_overlaps.Ranges <- function(x,y, maxgap = 0L, minoverlap =1L) {
  subsetByOverlaps(x,y,maxgap, minoverlap, type = "any",invert = TRUE)
}

filter_by_non_overlaps.GenomicRanges <- function(x,y, maxgap = 0L, minoverlap =1L) {
  subsetByOverlaps(x,y,maxgap, minoverlap, type = "any",
                   invert = TRUE,
                   ignore.strand = TRUE)
}
