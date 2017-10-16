#' Count the number of overlaps between two Ranges objects
#'
#' @param x,y Objects representing ranges
#' @param maxgap,minoverlap The maximimum gap between intervals as an integer
#' greater than or equal to zero. The minimum amount of overlap between intervals
#' as an integer greater than zero, accounting for the maximum gap.
#'
#' @return An integer vector of same length as x.
#' @importFrom IRanges countOverlaps
#' @rdname ranges-count-overlaps.Rd
#' @export
count_overlaps <- function(x, y, maxgap, minoverlap) {
  UseMethod("count_overlaps")
}

count_overlaps.Ranges <- function(x,y, maxgap = 0L, minoverlap = 1L) {
  countOverlaps(x,y, maxgap, minoverlap, type = "any")
}

count_overlaps.GenomicRanges <- function(x,y, maxgap = 0L, minoverlap = 1L) {
  countOverlaps(x,y, maxgap, minoverlap, type = "any", ignore.strand = TRUE)
}

#' @rdname ranges-count-overlaps.Rd
#' @export
count_overlaps_within <- function(x, y, maxgap, minoverlap) {
  UseMethod("count_overlaps_within")
}

count_overlaps_within.Ranges <- function(x,y, maxgap = 0L, minoverlap = 1L) {
  countOverlaps(x,y, maxgap, minoverlap, type = "within")
}

count_overlaps_within.GenomicRanges <- function(x,y, maxgap = 0L, minoverlap = 1L) {
  countOverlaps(x,y, maxgap, minoverlap, type = "within", ignore.strand = TRUE)
}
