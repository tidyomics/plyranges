# # ranges-anchors.R
#' Anchored Ranges objects
#'
#' @description The `GRangesAnchored` class and the `IRangesAnchored`
#' class allow components of a `GRanges` or `IRanges` (start, end, center)
#' to be held fixed.
#'
#' @details Anchoring will fix a Ranges start, end, or center positions,
#' so these positions will remain the same when performing arithimetic.
#' For `GRanges` objects, the function
#' (`anchor_3p()`) will fix the start for the negative strand,
#' while `anchor_5p()` will fix the end for the
#' positive strand. Anchoring modifies how arithmetic is performed, for example
#' modifying the width of a range with `set_width()` or stretching a
#' range with `stretch()`. To remove anchoring use `unanchor()`.
#'
#' @param x a Ranges object
#'
#' @section Constructors:
#' Depending on how you want to fix the components of a Ranges, there are
#' five ways to construct a RangesAnchored class. Here `x` is either
#' an `IRanges` or `GRanges` object.
#' \itemize{
#'    \item{`anchor_start(x)`}{Fix the start coordinates}
#'    \item{`anchor_end(x)`}{Fix the end coordinates}
#'    \item{`anchor_center(x)`}{Fix the center coordinates}
#'    \item{`anchor_3p(x)`}{On the negative strand fix the start coordinates,
#'    and for positive or unstranded ranges fix the end coordinates.}
#'    \item{`anchor_5p(x)`}{On the positive or unstranded ranges fix the start coordinates,
#'    coordinates and for negative stranded ranges fix the end coordinates.}
#' }
#'
#' @section Accessors:
#' To see what has been anchored use the function `anchor`.
#' This will return a character vector containing a valid anchor.
#' It will be set to one of `c("start", "end", "center")` for an
#' `IRanges` object or one of
#' `c("start", "end", "center", "3p", "5p")` for a `GRanges` object.
#'
#' @seealso \link{set_width}, \link{stretch}
#'
#' @return a RangesAnchored object which has the same appearance as a regular
#' Ranges object but with an additional slot displaying an anchor.
#' @examples
#' df <- data.frame(start = 1:10, width = 5)
#' rng <- as_iranges(df)
#' rng_by_start <- anchor_start(rng)
#' rng_by_start
#' anchor(rng_by_start)
#' set_width(rng_by_start, 3L)
#' grng <- as_granges(df,
#'                    seqnames = "chr1",
#'                    strand = c(rep("-", 5), rep("+", 5)))
#' rng_by_5p <- anchor_5p(grng)
#' rng_by_5p
#' set_width(rng_by_5p, 3L)
#'
#' @importFrom methods setClass setValidity setMethod show
#' @rdname ranges-anchor
#' @export
anchor <- function(x) { UseMethod("anchor") }

#' @export
anchor.AnchoredGenomicRanges <- function(x) {
  x@anchor
}

#' @export
anchor.AnchoredIntegerRanges <- anchor.AnchoredGenomicRanges

#' @rdname ranges-anchor
#' @export
unanchor <- function(x) { UseMethod("unanchor") }

#' @export
unanchor.AnchoredGenomicRanges <- function(x) {
  x@delegate
}

#' @export
unanchor.AnchoredIntegerRanges <- unanchor.AnchoredGenomicRanges


#' @rdname ranges-anchor
#' @export
anchor_start <- function(x) { UseMethod("anchor_start") }

#' @export
anchor_start.IntegerRanges <- function(x) {
  new_anchored_ir(x, "start")
}

#' @export
anchor_start.GenomicRanges <- function(x) {
  new_anchored_gr(x, "start")
}


#' @rdname ranges-anchor
#' @export
anchor_end <- function(x) { UseMethod("anchor_end") }

#' @export
anchor_end.IntegerRanges <- function(x) {
  new_anchored_ir(x, "end")
}

#' @export
anchor_end.GenomicRanges <- function(x) {
  new_anchored_gr(x, anchor = "end")
}

#' @rdname ranges-anchor
#' @export
anchor_center <- function(x) { UseMethod("anchor_center") }

#' @export
anchor_center.IntegerRanges <- function(x) {
  new_anchored_ir(x, "center")
}

#' @export
anchor_center.GenomicRanges <- function(x) {
  new_anchored_gr(x, "center")
}

#' @rdname ranges-anchor
#' @export
anchor_centre <- function(x) { UseMethod("anchor_center") }

#' @rdname ranges-anchor
#' @export
anchor_3p <- function(x) { UseMethod("anchor_3p") }

#' @export
anchor_3p.GenomicRanges <- function(x) {
  new_anchored_gr(x, "3p")
}

#' @rdname ranges-anchor
#' @export
anchor_5p <- function(x) { UseMethod("anchor_5p") }

#' @export
anchor_5p.GenomicRanges <- function(x) {
  new_anchored_gr(x, "5p")
}
