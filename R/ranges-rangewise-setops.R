# ranges-rangewise-setops.R

#' Vector-wise Range set-opterations
#'
#' @param x,y Two Ranges objects to compare.
#'
#' @details These are usual set-operations that act on the sets of the
#' ranges represented in x and y. By default these operations will ignore
#' any strand information. The directed versions of these functions will
#' take into account strand.
#' @return A Ranges object
#'
#' @examples
#' gr1 <- data.frame(seqnames = "chr1",
#'                   start = c(2,9),
#'                   end = c(7,9),
#'                   strand = c("+", "-")) %>%
#'                as_granges()
#' gr2 <- data.frame(seqnames = "chr1", start = 5, width = 5, strand = "-") %>%
#'          as_granges()
#'
#' union_ranges(gr1, gr2)
#' union_ranges_directed(gr1, gr2)
#'
#' intersect_ranges(gr1, gr2)
#' intersect_ranges_directed(gr1, gr2)
#'
#' setdiff_ranges(gr1, gr2)
#' setdiff_ranges_directed(gr1, gr2)
#'
#' @rdname ranges-setops
#' @export
intersect_ranges <- function(x,y) { UseMethod("intersect_ranges") }

#' @importFrom IRanges intersect
#' @export
intersect_ranges.IntegerRanges <- function(x,y) {
  intersect(x,y)
}

#' @importFrom GenomicRanges intersect
#' @export
intersect_ranges.GenomicRanges <- function(x,y) {
  intersect(x,y,ignore.strand = TRUE)
}

#' @export
#' @rdname ranges-setops
intersect_ranges_directed <- function(x,y) { UseMethod("intersect_ranges_directed") }

#' @export
#' @importFrom GenomicRanges intersect
intersect_ranges_directed.GenomicRanges <- function(x,y) {
  intersect(x,y, ignore.strand = FALSE)
}

#' @export
#' @rdname ranges-setops
union_ranges <- function(x,y) { UseMethod("union_ranges") }

#' @export
#' @importFrom IRanges union
union_ranges.IntegerRanges <- function(x,y) {
  union(x,y)
}

#' @export
#' @importFrom GenomicRanges union
union_ranges.GenomicRanges <- function(x,y) {
  union(x,y,ignore.strand = TRUE)
}

#' @export
#' @rdname ranges-setops
union_ranges_directed <- function(x,y) { UseMethod("union_ranges_directed") }

#' @export
#' @importFrom GenomicRanges union
union_ranges_directed.GenomicRanges <- function(x,y) {
  union(x,y, ignore.strand = FALSE)
}

#' @export
#' @rdname ranges-setops
setdiff_ranges <- function(x,y) { UseMethod("setdiff_ranges") }

#' @export
#' @importFrom IRanges setdiff
setdiff_ranges.IntegerRanges <- function(x,y) {
  setdiff(x,y)
}

#' @export
#' @importFrom GenomicRanges setdiff
setdiff_ranges.GenomicRanges <- function(x,y) {
  setdiff(x,y, ignore.strand = TRUE)
}

#' @export
#' @rdname ranges-setops
setdiff_ranges_directed <- function(x,y) { UseMethod("setdiff_ranges_directed") }

#' @export
#' @importFrom GenomicRanges setdiff
setdiff_ranges_directed.GenomicRanges <- function(x,y) {
  setdiff(x,y, ignore.strand = FALSE)
}
