# ranges-rangewise-setops.R

#' Vector-wise Range set-opterations
#'
#' @param x,y Two Ranges objects to compare.
#'
#' @details These are usual set-operations that act on the sets of the
#' ranges represented in x and y. By default these operations will ignore
#' any strand information. The directedd versions of these functions will
#' take into account strand.
#'
#' @export
#' @rdname ranges-setops
intersect_ranges <- function(x,y) { UseMethod("intersect_ranges") }

#' @importFrom IRanges intersect
intersect_ranges.Ranges <- function(x,y) {
  intersect(x,y)
}

#' @importFrom GenomicRanges intersect
intersect_ranges.GenomicRanges <- function(x,y) {
  intersect(x,y,ignore.strand = TRUE)
}

#' @export
#' @rdname ranges-setops
intersect_ranges_directed <- function(x,y) { UseMethod("intersect_ranges_directed") }

#' @importFrom GenomicRanges intersect
intersect_ranges_directed.GenomicRanges <- function(x,y) {
  intersect(x,y, ignore.strand = FALSE)
}

#' @export
#' @rdname ranges-setops
union_ranges <- function(x,y) { UseMethod("union_ranges") }

#' @importFrom IRanges union
union_ranges.Ranges <- function(x,y) {
  union(x,y)
}

#' @importFrom GenomicRanges union
union_ranges.GenomicRanges <- function(x,y) {
  union(x,y,ignore.strand = TRUE)
}

#' @export
#' @rdname ranges-setops
union_ranges_directed <- function(x,y) { UseMethod("union_ranges_directed") }

#' @importFrom GenomicRanges union
union_ranges_directed.GenomicRanges <- function(x,y) {
  union(x,y, ignore.strand = FALSE)
}

#' @export
#' @rdname ranges-setops
setdiff_ranges <- function(x,y) { UseMethod("setdiff_ranges") }

#' @importFrom IRanges setdiff
setdiff_ranges.Ranges <- function(x,y) {
  setdiff(x,y)
}

#' @importFrom GenomicRanges setdiff
setdiff_ranges.Ranges <- function(x,y) {
  setdiff(x,y, ignore.strand = FALSE)
}

#' @export
#' @rdname ranges-setops
setdiff_ranges_directed <- function(x,y) { UseMethod("setdiff_ranges_directed") }

#' @importFrom GenomicRanges setdiff
setdiff_ranges_directed.GenomicRanges <- function(x,y) {
  setdiff(x,y, ignore.strand = FALSE)
}
