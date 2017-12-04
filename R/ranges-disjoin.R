# ranges-disjoin.R
#' Disjoin then aggregate a Ranges object
#'
#' @param .data a Ranges object to disjoin
#' @param ... Name-value pairs of summary functions.
#'
#' @return a Ranges object with the
#' @rdname ranges-disjoin
#' @importFrom IRanges disjoin
#' @examples
#' df <- data.frame(start = 1:10, width = 5,  seqnames = "seq1",
#' strand = sample(c("+", "-", "*"), 10, replace = TRUE), gc = runif(10))
#' rng <- Ranges(df)
#' rng %>% disjoin_ranges()
#' rng %>% disjoin_ranges(gc = mean(gc))
#' rng %>% disjoin_ranges_directed(gc = mean(gc))
#' @export
disjoin_ranges <- function(.data, ...) { UseMethod("disjoin_ranges") }

#' @method disjoin_ranges Ranges
#' @export
disjoin_ranges.Ranges <- function(.data, ...) {
  dots <- quos(...)
  if (length(dots) == 0L) {
    return(disjoin(.data))
  }

  disjoined <- disjoin(.data, with.revmap = TRUE)

  reduce_rng(.data, disjoined, dots)

}

#' @method disjoin_ranges GenomicRanges
#' @export
disjoin_ranges.GenomicRanges <- function(.data, ...) {

  dots <- quos(...)
  if (length(dots) == 0L) {
    return(disjoin(.data,ignore.strand = TRUE))
  }

  disjoined <- disjoin(.data, with.revmap = TRUE, ignore.strand = TRUE)
  reduce_rng(.data, disjoined, dots)
}

#' @rdname ranges-disjoin
#' @export
disjoin_ranges_directed <- function(.data, ...) {
  UseMethod("disjoin_ranges_directed")
}

#' @importFrom IRanges disjoin
#' @method disjoin_ranges_directed GenomicRanges
#' @export
disjoin_ranges_directed.GenomicRanges <- function(.data, ...) {

  dots <- quos(...)
  if (length(dots) == 0L) {
    return(disjoin(.data,ignore.strand = FALSE))
  }

  disjoined <- disjoin(.data, with.revmap = TRUE, ignore.strand = FALSE)
  reduce_rng(.data, disjoined, dots)

}

