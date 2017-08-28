# between-joins.R
# join verbs for Ranges-like objects without overlapping.
# These are wrapper functions to `merge` like GenomicRanges overlaps methods,
# and provide a consistent API for performing joins between the following
# classes:
# Ranges, DataFrame
# should there by a Ranges, tibble join?

#' Join two GRanges/tbls together
#' @name join
#' @rdname join
#' @description These are generic methods to perform basic joins
#' on combinations of GRanges objects and/or table-like objects.
#' For GRanges objects these joins do not take into account
#' genomic intervals but merely join the mcols tbls in a sensible way.
#' @param x,y
#' @param by a character vector of variables to join. The default is to perform
#' a natural join. Use a named vector to join by different variables on x and y.
#'
#' @seealso \link[dplyr]{join}
#' @importFrom S4Vectors merge
#' @importFrom BiocGenerics intersect
#' @importFrom GenomicRanges mcols
#' @export
setGeneric("left_join",
           function(x, y, by = NULL, ...) standardGeneric("left_join"))
#' @describeIn join
#' @export
setGeneric("right_join",
           function(x, y, by = NULL, ...) standardGeneric("right_join"))
#' @describeIn join
#' @export
setGeneric("inner_join",
           function(x, y, by = NULL, ...) standardGeneric("inner_join"))
#' @describeIn join
#' @export
setGeneric("full_join",
           function(x, y, by = NULL, ...) standardGeneric("full_join"))

#' @export
setGeneric("common_mcols",
           function(x, y, by = NULL) standardGeneric("common_mcols"))

# should perform type checking, but does not
#' @export
setMethod("common_mcols", c("DataFrame", "DataFrame"),
          function(x, y, by = NULL) {
            x_names <- names(x)
            y_names <- names(y)
            if (is.null(by)) {
              common_mcols <- intersect(x_names,  y_names)
              if (length(common_mcols) == 0) {
                stop("No common columns between x & y", call. = FALSE)
              }
              return(common_mcols)
            } else {
              named_by <- names(by)
              if (length(named_by) > 0) {
                stopifnot(named_by %in% x_names || by %in% y_names)
                by

              } else {
                stopifnot(by %in% x_names || by %in% y_names)
                by
              }
            }

          })

#' @export
setMethod("common_mcols", c("GRanges", "GRanges"),
          function(x, y, by = NULL) {
            common_mcols(mcols(x), mcols(y), by = by)
          })

#' @describeIn join
#' @export
setMethod("inner_join", c("DataFrame", "DataFrame"),
          function(x, y, by = NULL) {

            join_by <- common_mcols(x,y, by)
            is_named_by <- length(names(join_by)) > 0

            if (is_named_by) {
              y_col <- join_by
              x_col <- names(join_by)
              merge(x, y, by.x = x_col, by.y = y_col,
                    by = NULL, sort = FALSE)
            } else {
              merge(x, y, by = join_by, sort = FALSE)
            }

          })

# an inner join reduces the metadata DataFrames to their common rows
# if we naively copy that to the ranges in x we will get an error,
# so we need to subset the ranges in x based on the keys in the by
# Note this is different to an overlap join - we are just joining the
# mcols(x), mcols(y)
#' @describeIn join
#' @export
setMethod("inner_join", c("GRanges", "GRanges"),
          function(x, y, by = NULL) {
            join_by <- common_mcols(x,y, by = by)
            is_named <- length(names(join_by)) > 0
            if (!is_named) {
              names(join_by) <- join_by
            }

            ranges_common <- x[filter_common(join_by, x, y)]
            mcols(ranges_common) <- inner_join(mcols(ranges_common),
                                               mcols(y), by = by)
            ranges_common
          })

#' @describeIn join
#' @export
setMethod("inner_join", c("GRanges", "DataFrame"),
          function(x, y, by = NULL) {
            join_by <- common_mcols(mcols(x), y, by = by)
            is_named <- length(names(join_by)) > 0
            if (!is_named) {
              names(join_by) <- join_by
            }
            ranges_common <- x[filter_common(join_by, x, y)]
            mcols(ranges_common) <- inner_join(mcols(x), y)
            ranges_common
          })

#' @describeIn join
#' @export
setMethod("inner_join", c("DataFrame", "GRanges"),
          function(x,y, by = NULL) {
            inner_join(x, mcols(y), by = by)
          })

#' @describeIn join
#' @export
setMethod("left_join", c("DataFrame", "DataFrame"),
          function(x, y, by = NULL) {
            join_by <- common_mcols(x,y, by)
            is_named_by <- length(names(join_by)) > 0

            if (is_named_by) {
              y_col <- join_by
              x_col <- names(join_by)
              merge(x, y, by.x = x_col, by.y = y_col,
                    by = NULL, all.x = TRUE, sort = FALSE)
            } else {
              merge(x, y, by = join_by, all.x = TRUE, sort = FALSE)
            }
          })

# should have an option to copy here?
#' @describeIn join
#' @export
setMethod("left_join", c("GRanges", "GRanges"),
          function(x, y, by = NULL) {
            mcols(x) <- left_join(mcols(x), mcols(y), by = by)
            x
          })

#' @describeIn join
#' @export
setMethod("left_join", c("GRanges", "DataFrame"),
          function(x, y, by = NULL) {
            mcols(x) <- left_join(mcols(x), y, by = by)
            x
          })

#' @describeIn join
#' @export
setMethod("left_join", c("DataFrame", "GRanges"),
          function(x, y, by = NULL) {
            left_join(x, mcols(y), by = by)
          })
#' @describeIn join
#' @export
setMethod("right_join", c("DataFrame", "DataFrame"),
          function(x, y, by = NULL) {
            join_by <- common_mcols(x,y, by)
            is_named_by <- length(names(join_by)) > 0

            if (is_named_by) {
              y_col <- join_by
              x_col <- names(join_by)
              merge(x, y, by.x = x_col, by.y = y_col,
                    by = NULL, all.y = TRUE, sort = FALSE)
            } else {
              merge(x, y, by = join_by, all.y = TRUE, sort = FALSE)
            }

          })

#' @describeIn join
#' @export
setMethod("right_join", c("GRanges", "GRanges"),
          function(x, y, by = NULL) {
            mcols(y) <- right_join(mcols(x), mcols(y), by = by)
            y
          })

#' @describeIn join
#' @export
setMethod("right_join", c("GRanges", "DataFrame"),
          function(x, y, by = NULL) {
            right_join(mcols(x), y, by = by)
          })

#' @describeIn join
#' @export
setMethod("right_join", c("DataFrame", "GRanges"),
          function(x, y, by = NULL) {
            mcols(y) <- right_join(x, mcols(y), by = by)
            y
          })

#' @describeIn join
#' @export
setMethod("full_join", c("DataFrame", "DataFrame"),
          function(x, y, by = NULL) {
            join_by <- common_mcols(x,y, by)
            is_named_by <- length(names(join_by)) > 0

            if (is_named_by) {
              y_col <- join_by
              x_col <- names(join_by)
              merge(x, y, by.x = x_col, by.y = y_col,
                    by = NULL, all = TRUE, sort = FALSE)
            } else {
              merge(x, y, by = join_by, all = TRUE, sort = FALSE)
            }
          })

# unclear how the full join should work in terms of returning the ranges
# at the moment full join performs a full join on the metadata
# then returns the ranges on either x or y
#' @describeIn join
#' @export
setMethod("full_join", c("GRanges", "GRanges"),
          function(x,y, by = NULL) {

            join_by <- common_mcols(x,y, by = by)
            is_named <- length(names(join_by)) > 0
            if (!is_named) {
              names(join_by) <- join_by
            }

            rows_filter_y <- negate_rows_common(join_by, y, x)
            reduced_y <- y[Reduce(all, rows_filter_y)]
            mcols(reduced_y) <- NULL
            full_ranges <- c(x, reduced_y, ignore.mcols = TRUE)
            mcols(full_ranges) <- full_join(mcols(x), mcols(y), by = by)
            full_ranges

          })

#' @describeIn join
#' @export
setMethod("full_join", c("DataFrame", "GRanges"),
          function(x, y, by = NULL) {
            full_join(x, mcols(y), by = by)
          })
# does it make sense to return a Ranges here
#' @describeIn join
#' @export
setMethod("full_join", c("GRanges", "DataFrame"),
          function(x, y, by = NULL) {
            full_join(mcols(x), y, by = by)
          })
