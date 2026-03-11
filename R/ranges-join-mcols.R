.join_mcols <- function(x, y, by = NULL, all) {
  x_id <- x
  x_id$.id <- factor(seq_along(x))
  if (!is(y, "DataFrame")) {
    y <- DataFrame(y)
  }
  if (is.null(by)) {
    by <- intersect(names(mcols(x_id)), names(y))
    if (length(by) == 0) {
      rlang::abort(
        "`by` must be supplied when `x` and `y` have no common variables."
      )
    }
    message("Joining with `by = ", deparse(by), "`")
  }
  new_mcols <- S4Vectors::merge(mcols(x_id), y, by = by, sort = FALSE, all = all)
  new_mcols <- new_mcols[order(new_mcols$.id), ]
  new_x <- x[as.integer(new_mcols$.id)]
  new_mcols$.id <- NULL
  mcols(new_x) <- NULL
  mcols(new_x) <- new_mcols
  return(new_x)
}


#' Join data by metadata columns
#'
#' These functions enable joins based on metadata
#' of `x` and data in `y`, utilizing
#' [S4Vectors::merge()].
#'
#' @param x Object representing ranges, with metadata
#' columns containing variables for matching
#' @param y A table of data, DataFrame, data.frame or tibble
#' @param by Specifications of the columns used for merging.
#' Passed to `S4Vectors::merge()`
#' @param ... Additional arguments passed to `S4Vectors::merge()`
#' 
#' @details The function [join_mcols_inner()] returns ranges in 
#' `x` only when a match is found in the columns specified with
#' `by` in the table `y`.
#' 
#' The function [join_mcols_left()] returns all ranges in `x`
#' regardless of a match in `y`, with duplications possible
#' from multiple matches.
#' 
#' @return Object representing ranges, with new metadata columns
#'
#' @examples
#'
#' x <- as_granges(data.frame(
#'   seqnames = "chr1", start = 1:5, width=10,
#'   id=c(1:4,4), foo=letters[1:5]
#' ))
#' y <- DataFrame(id=c(1:2,2,4), bar=letters[6:9])
#'
#' # metadata joins:
#' join_mcols_inner(x, y, by="id")
#' join_mcols_left(x, y, by="id")
#'
#' @importFrom dplyr left_join
#' @importFrom rlang .data
#' @importFrom S4Vectors merge
#'
#' @rdname mcols-joins
#' @export
join_mcols_inner <- function(x, y, by = NULL, ...) {
  .join_mcols(x = x, y = y, by = by, all = FALSE, ...)
}

#' @rdname mcols-joins
#' @export
join_mcols_left <- function(x, y, by = NULL, ...) {
  .join_mcols(x = x, y = y, by = by, all = TRUE, ...)
}
