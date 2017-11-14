#' Find overlaps within a Ranges object
#'
#' @param x A Ranges object
#' @param maxgap,minoverlap The maximimum gap between intervals as an integer
#' greater than or equal to zero. The minimum amount of overlap between intervals
#' as an integer greater than zero, accounting for the maximum gap.
#'
#' @details Self overlaps find any overlaps (or overlaps contained) within a
#' ranges and append the overlapping details as metadata columns. A self inner
#' overlap join will update the ranges object with all self intersecting ranges.
#'
#' @rdname self-overlap-joins
#' @export
join_overlap_self <- function(x, maxgap, minoverlap) {
  UseMethod("join_overlap_self")
}

join_overlap_self.Ranges <- function(x, maxgap = -1L, minoverlap = 0L) {
  find_overlaps(x,x, maxgap, minoverlap, suffix = c("", ".overlap"))
}

join_overlap_self.GenomicRanges <- function(x, maxgap = -1L, minoverlap = 0L) {
  find_overlaps(x,x, maxgap, minoverlap, suffix = c("", ".overlap"))
}

#' @rdname self-overlap-joins
#' @export
join_overlap_within_self <- function(x, maxgap, minoverlap) {
  UseMethod("join_overlap_within_self")
}

join_overlap_within_self.Ranges <- function(x, maxgap = -1L, minoverlap = 0L) {
  find_overlaps_within(x,x, maxgap, minoverlap, suffix = c("", ".overlap"))
}

join_overlap_within_self.GenomicRanges <- function(x, maxgap = -1L, minoverlap = 0L) {
  find_overlaps_within(x,x, maxgap, minoverlap, suffix = c("", ".overlap"))
}


#' @rdname self-overlap-joins
#' @export
join_overlap_self_inner <- function(x, maxgap, minoverlap) {
  UseMethod("join_overlap_within_self")
}

join_overlap_self_inner.Ranges <- function(x, maxgap = 0L, minoverlap = 1L) {
  join_overlap_inner(x,x, maxgap, minoverlap, suffix = c("", ".overlap"))
}

join_overlap_within_self.GenomicRanges <- function(x, maxgap = 0L, minoverlap = 1L) {
  join_overlap_inner(x,x, maxgap, minoverlap, suffix = c("", ".overlap"))
}
