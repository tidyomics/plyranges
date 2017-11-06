# less error prone to just import dplyr's joins here?
join_left <- function(x, y, by = NULL) {
  join_by <- common_cols(x,y, by)
  is_named_by <- length(names(join_by)) > 0

  x_os <- overscope_ranges(x)
  on.exit(overscope_clean(x_os))
  y_os <- overscope_ranges(y)
  on.exit(overscope_clean(y_os), add = TRUE)

  if (is_named_by) {
    y_col <- syms(join_by)
    names(y_col) <- join_by
    x_col <- syms(names(join_by))
    names(x_col) <- names(join_by)

  } else {
    x_col <- y_col <- syms(join_by)
    names(x_col) <- names(y_col) <- join_by
  }

  x_list <- lapply(x_col, overscope_eval_next,  overscope = x_os)
  y_list <- lapply(y_col, overscope_eval_next, overscope = y_os)
  inx <- lapply(seq_along(join_by), function(i) {
    x_col <- as.character(x_col[[i]])
    y_col <- as.character(y_col[[i]])
    x <- BiocGenerics::match(x_list[[x_col]], y_list[[y_col]])
    x[!is.na(x)]
  })
  inx <- Reduce(intersect, inx)
  if (length(inx) == 0) {
    stop("No common keys between x & y", .call = FALSE)
  }
  # suffix here
  mcols(x) <- cbind(mcols(x), mcols(y[inx]))
  return(x)
}

join_inner <- function(x, y, by = NULL) {
  join_by <- common_cols(x,y, by)
  is_named_by <- length(names(join_by)) > 0

  x_os <- overscope_ranges(x)
  on.exit(overscope_clean(x_os))
  y_os <- overscope_ranges(y)
  on.exit(overscope_clean(y_os), add = TRUE)

  if (is_named_by) {
    y_col <- syms(join_by)
    names(y_col) <- join_by
    x_col <- syms(names(join_by))
    names(x_col) <- names(join_by)

  } else {
    x_col <- y_col <- syms(join_by)
    names(x_col) <- names(y_col) <- join_by
  }

  x_list <- lapply(x_col, overscope_eval_next,  overscope = x_os)
  y_list <- lapply(y_col, overscope_eval_next, overscope = y_os)
  inx <- lapply(seq_along(join_by), function(i) {
    x_col <- as.character(x_col[[i]])
    y_col <- as.character(y_col[[i]])
    list(lhs = BiocGenerics::match(x_list[[x_col]], y_list[[y_col]]),
         rhs = BiocGenerics::match(y_list[[y_col]], x_list[[x_col]]))
  })

  inx_x <- Reduce(intersect, lapply(inx, function(x) x$lhs[!is.na(x$lhs)]))
  inx_y <- Reduce(intersect, lapply(inx, function(x) x$rhs[!is.na(x$rhs)]))

  ranges_common <- x[inx_x][inx_y]
  mcols(ranges_common) <- cbind(mcols(ranges_common), mcols(y[seq_along(inx_y)]))
  ranges_common

}


#' Join metadata from Ranges together
#' @name join
#' @rdname join
#' @description These are methods to perform basic joins
#' on combinations of Ranges objects and/or table-like objects.
#' For Ranges objects these joins do not take into account
#' genomic intervals but merely join the metadata columns in a sensible way.
#' @param x,y Ranges objects to join together
#' @param by a character vector of variables to join. The default is to perform
#' a natural join. Use a named vector to join by different variables on x and y.
#'
#' @seealso \link[dplyr]{join}
#' @importFrom BiocGenerics intersect
#' @importFrom GenomicRanges mcols
#' @importFrom S4Vectors merge
#' @importFrom dplyr inner_join
#' @method inner_join Ranges
#' @examples
#' set.seed(100)
#' df_a <- data.frame(start = 1:10, width = 5,
#'                  id = rep(letters[1:2], 5),
#'                  id2 = sample(letters[1:2], 10, replace = TRUE))
#' df_b <- data.frame(end = 5:14, width = 10,
#'                    id = sort(rep(letters[1:2], 5)),
#'                    z = rnorm(10))
#' rng_a <- Ranges(df_a)
#' rng_b <- Ranges(df_b)
#' inner_join(rng_a, rng_b, by = "id")
#' inner_join(rng_a, rng_b, by = c("id2" = "id"))
#' left_join(rng_a, rng_b, by = "id")
#' left_join(rng_b, rng_a, by = "id")
#' @export
inner_join.Ranges <- function(x, y, by = NULL) {
  stopifnot(is_ranges(y))
  join_inner(x, y, by)
}

#' @method inner_join GenomicRanges
#' @export
inner_join.GenomicRanges <- function(x, y, by = NULL) {
  stopifnot(is_ranges(y))
  join_inner(x, y, by)
}


#' @importFrom S4Vectors merge
#' @importFrom dplyr left_join
#' @export
left_join.Ranges <- function(x, y, by = NULL) {
  stopifnot(is_ranges(y))
  join_left(x, y, by)
}

#' @method left_join GenomicRanges
#' @export
left_join.GenomicRanges <- function(x, y, by = NULL) {
  stopifnot(is_ranges(y))
  join_left(x, y, by)
}
