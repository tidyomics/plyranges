
#' Filter by overlapping/non-overlapping ranges
#'
#' @param x,y Objects representing ranges
#' @param maxgap The maximimum gap between intervals as a single
#' integer greater than or equal to -1. If you modify this argument,
#' \code{minoverlap} must be held fixed.
#' @param minoverlap  The minimum amount of overlap between intervals
#' as a single integer greater than 0. If you modify this argument,
#' \code{maxgap} must be held fixed.
#'
#' @details By default, \code{filter_by_overlaps} and
#' \code{filter_by_non_overlaps} ignore strandedness for \code{\link{GRanges}}
#' objects. The argument \code{maxgap} is the maximum number of positions
#' between two ranges for them to be considered overlapping. Here the default
#' is set to be -1 as that is the the gap between two ranges that
#' has its start or end strictly inside the other. The argugment
#' \code{minoverlap} refers to the minimum number of positions
#' overlapping between ranges, to consider there to be overlap.
#'
#' @importFrom IRanges subsetByOverlaps
#' @seealso \link[IRanges]{subsetByOverlaps}
#' @export
#' @examples
#' df <- data.frame(seqnames = c("chr1", rep("chr2", 2),
#'                               rep("chr3", 3), rep("chr4", 4)),
#'                  start = 1:10,
#'                  width = 10:1,
#'                  strand = c("-", "+", "+", "*", "*", "+", "+", "+", "-", "-"),
#'                  name = letters[1:10])
#' query <- as_granges(df)
#'
#' df2 <- data.frame(seqnames = c(rep("chr2", 2), rep("chr1", 3), "chr2"),
#'                   start = c(4,3,7,13,1,4),
#'                   width = c(6,6,3,3,3,9),
#'                   strand = c(rep("+", 3), rep("-", 3)))
#' subject <- as_granges(df2)
#'
#' filter_by_overlaps(query, subject)
#'
#' filter_by_non_overlaps(query, subject)
#'
#' @rdname filter-overlaps.rd
filter_by_overlaps <- function(x,y, maxgap = -1L, minoverlap = 0L) {
  UseMethod("filter_by_overlaps")
}

#' @export
filter_by_overlaps.Ranges <- function(x, y, maxgap = -1L, minoverlap = 0L) {
  subsetByOverlaps(x,y, maxgap, minoverlap)
}

#' @export
filter_by_overlaps.GenomicRanges <- function(x, y, maxgap = -1L, minoverlap =0L) {
  subsetByOverlaps(x,y, maxgap, minoverlap, ignore.strand = TRUE)
}

#' @export
#' @rdname filter-overlaps.rd
filter_by_non_overlaps <- function(x,y, maxgap, minoverlap) {
  UseMethod("filter_by_non_overlaps")
}

#' @export
filter_by_non_overlaps.Ranges <- function(x,y, maxgap = -1L, minoverlap = 0L) {
  subsetByOverlaps(x,y, maxgap, minoverlap, invert = TRUE)
}

#' @export
filter_by_non_overlaps.GenomicRanges <- function(x,y, maxgap = -1L, minoverlap = 0L) {
  subsetByOverlaps(x,y, maxgap, minoverlap,
                   invert = TRUE,
                   ignore.strand = TRUE)
}
