# ranges-join-follow.R

#' Find following Ranges
#'
#' @param x,y Ranges objects, which ranges in x follow those in y.
#'
#' @details By default \code{join_follow} will find abritrary ranges
#' in y that are followed by ranges in x and ignore any strand information.
#' On the other hand \code{join_follow_left} will find all ranges in y
#' that are on the left-hand side of the ranges in x ignoring any strand
#' information. Finally, \code{join_follow_upstream} will find all ranges in x
#' that are that are upstream of the ranges in y. On the positive strand this
#' will result in ranges in y that are left of those in x and on the negative
#' strand it will result in ranges in y that are right of those in x.
#'
#' @return A Ranges object with a metadata column called follows that
#' contains the corresponding Ranges in y that are followed by the ranges in x.
#'
#' @rdname follow-ranges
#' @importFrom IRanges follow
#' @export
join_follow <- function(x,y) { UseMethod("join_follow") }

#' @export
join_follow.Ranges <- function(x,y) {
  hits <- follow(x,y)
  nearest_rng(x,y,hits, type = "follows")
}

#' @export
join_follow.GenomicRanges <- function(x,y) {
  hits <- follow(x,y, ignore.strand = TRUE)
  nearest_rng(x,y, hits, type = "follows")
}



#' @rdname follow-ranges
#' @importFrom IRanges follow
#' @export
join_follow_left <- function(x,y) { UseMethod("join_follow_left") }

#' @export
join_follow_left.Ranges <- function(x,y) {
  hits <- follow(x,y, select = "all")
  nearest_rng_all(x,y, hits, type = "follows")
}

#' @export
join_follow_left.GenomicRanges <- function(x,y) {
  hits <- follow(x,y, select = "all", ignore.strand = TRUE)
  nearest_rng_all(x,y, hits, type = "follows")
}

#' @rdname follow-ranges
#' @importFrom IRanges follow
#' @export
join_follow_upstream <- function(x,y) {UseMethod("join_follow_upstream")}

#' @export
join_follow_upstream.GenomicRanges <- function(x,y) {
  hits <- follow(x,y, select = "all", ignore.strand = FALSE)
  nearest_rng_all(x, y, hits, type = "follows")
}
