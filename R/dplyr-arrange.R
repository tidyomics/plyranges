rng_arrange <- function(.data, dots) {
  rng_os <- overscope_ranges(.data)
  on.exit(overscope_clean(rng_os))
  rng_list <- lapply(dots, overscope_eval_next, overscope = rng_os)

  if (length(rng_list) == 1L) {
    idx <- order(rng_list[[1]])
  } else {
    idx <- Reduce(order, rng_list)
  }
  .data[idx, ]
}
#' Sort a Ranges object
#'
#' @param .data A Ranges object.
#' @param ... Comma seperated list of variable names.
#'
#' @importFrom dplyr arrange
#' @importFrom rlang quos
#' @examples
#' rng <- as_iranges(data.frame(start = 1:10, width = 10:1))
#' rng <- mutate(rng, score = runif(10))
#' arrange(rng, score)
#' # you can also use dplyr::desc to arrange by descending order
#' @rdname ranges-arrange
#' @return A sorted Ranges object
#' @method arrange Ranges
#' @export
arrange.Ranges <- function(.data, ...) {
  dots <- quos(...)
  rng_arrange(.data, dots)
}

#' @method arrange DelegatingGenomicRanges
#' @export
arrange.DelegatingGenomicRanges <- function(.data, ...) {
  dots <- quos(...)
  delegate <- .data@delegate
  .data@delegate <- rng_arrange(delegate, dots)
  .data
}

#' @method arrange DelegatingGenomicRanges
#' @export
arrange.DelegatingIntegerRanges <- function(.data, ...) {
  dots <- quos(...)
  delegate <- .data@delegate
  .data@delegate <- rng_arrange(delegate, dots)
  .data
}
# -- TODO: allow sorting by group?