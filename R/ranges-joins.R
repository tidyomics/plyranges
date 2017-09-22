# less error prone to just import dplyr's joins here?
join_left <- function(x, y, by = NULL) {
  join_by <- common_mcols(x,y, by)
  is_named_by <- length(names(join_by)) > 0

  if (is_named_by) {
    y_col <- join_by
    x_col <- names(join_by)
    mcols(x) <- merge(mcols(x), mcols(y), by.x = x_col, by.y = y_col,
                      by = NULL, all.x = TRUE, sort = FALSE)
  } else {
    mcols(x) <- merge(mcols(x), mcols(y), by = join_by, all.x = TRUE, sort = FALSE)
  }
  return(x)
}

join_inner <- function(x, y, by = NULL) {
  join_by <- common_mcols(x,y, by)
  is_named_by <- length(names(join_by)) > 0

  ranges_common <- x[filter_common(join_by, x, y)]

  if (is_named_by) {
    y_col <- join_by
    x_col <- names(join_by)
    mcols(ranges_common) <- merge(mcols(x),
                                  mcols(y),
                                  by.x = x_col, by.y = y_col,
                                  by = NULL, sort = FALSE)
  } else {
    mcols(ranges_common) <- merge(mcols(x), mcols(y), by = join_by, sort = FALSE)
  }
  return(ranges_common)
}


#' Join metadata from Ranges together
#' @name join
#' @rdname join
#' @description These are methods to perform basic joins
#' on combinations of Ranges objects and/or table-like objects.
#' For Ranges objects these joins do not take into account
#' genomic intervals but merely join the metadata columns in a sensible way.
#' @param x,y
#' @param by a character vector of variables to join. The default is to perform
#' a natural join. Use a named vector to join by different variables on x and y.
#'
#' @seealso \link[dplyr]{join}
#' @importFrom BiocGenerics intersect
#' @importFrom GenomicRanges mcols
#' @importFrom S4Vectors merge
#' @importFrom dplyr inner_join
#' @export
inner_join.Ranges <- function(x, y, by = NULL) {
  stopifnot(is_ranges(y))
  join_inner(x, y, by)
}

#' @export
inner_join.GenomicRanges <- function(x, y, by = NULL) {
  stopifnot(is_ranges(y))
  join_inner(x, y, by)
}


#' @importFrom S4Vectors merge
#' @importFrom dplyr inner_join
#' @export
left_join.Ranges <- function(x, y, by = NULL) {
  stopifnot(is_ranges(y))
  join_left(x, y, by)
}

#' @export
left_join.GenomicRanges <- function(x, y, by = NULL) {
  stopifnot(is_ranges(y))
  join_left(x, y, by)
}
