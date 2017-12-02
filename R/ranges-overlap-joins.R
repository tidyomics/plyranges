#' @rdname overlap-joins
#' @export
join_overlap_inner <- function(x, y, maxgap, minoverlap, suffix = c(".x", ".y")) {
  UseMethod("join_overlap_inner")
}

#' @rdname overlap-joins
#' @export
join_overlap_inner.Ranges <- function(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  find_overlaps(x,y,maxgap, minoverlap,suffix)
}

#' @rdname overlap-joins
#' @export
join_overlap_inner.GenomicRanges <- function(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  find_overlaps(x,y,maxgap, minoverlap,suffix)
}

#' @rdname overlap-joins
#' @export
join_overlap_inner_within <- function(x, y, maxgap, minoverlap, suffix = c(".x", ".y")) {
  UseMethod("join_overlap_inner_within")
}

#' @rdname overlap-joins
#' @export
join_overlap_inner_within.Ranges <- function(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  find_overlaps_within(x,y,maxgap, minoverlap,suffix)
}

#' @rdname overlap-joins
#' @export
join_overlap_inner_within.GenomicRanges <- function(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  find_overlaps_within(x,y, maxgap, minoverlap,suffix)
}

#' @rdname overlap-joins
#' @export
join_overlap_inner_directed <- function(x, y, maxgap, minoverlap, suffix = c(".x", ".y")) {
  UseMethod("join_overlap_inner_directed")
}

#' @rdname overlap-joins
#' @export
join_overlap_inner_directed.GenomicRanges <- function(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  find_overlaps_directed(x,y, maxgap, minoverlap,suffix)
}
