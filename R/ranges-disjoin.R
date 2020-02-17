#' Disjoin then aggregate a Ranges object
#'
#' @param .data a Ranges object to disjoin
#' @param ... Name-value pairs of summary functions.
#'
#' @return a Ranges object that is now disjoint (no bases overlap). 
#' @rdname ranges-disjoin
#' @importFrom IRanges disjoin
#' @examples
#' df <- data.frame(start = 1:10, width = 5,  seqnames = "seq1",
#' strand = sample(c("+", "-", "*"), 10, replace = TRUE), gc = runif(10))
#' rng <- as_granges(df)
#' rng %>% disjoin_ranges()
#' rng %>% disjoin_ranges(gc = mean(gc))
#' rng %>% disjoin_ranges_directed(gc = mean(gc))
#' @export
disjoin_ranges <- function(.data, ...) { UseMethod("disjoin_ranges") }

#' @method disjoin_ranges IntegerRanges
#' @export
disjoin_ranges.IntegerRanges <- function(.data, ...) {
  reduce_single(.data, ..., rfun = disjoin)
}

#' @method disjoin_ranges GroupedIntegerRanges
#' @export
disjoin_ranges.GroupedIntegerRanges <- function(.data, ...) {
  reduce_by_grp(.data, ..., rfun = disjoin)
}

#' @method disjoin_ranges GenomicRanges
#' @export
disjoin_ranges.GenomicRanges <- function(.data, ...) {
  reduce_single(.data, ...,
                rfun = function(x, ...) {
                  disjoin(x, ..., ignore.strand = TRUE)
                })
}

#' @method disjoin_ranges GroupedGenomicRanges
#' @export
disjoin_ranges.GroupedGenomicRanges <- function(.data, ...) {
  reduce_by_grp(.data, ...,
                rfun = function(x, ...) {
                  disjoin(x, ..., ignore.strand = TRUE)
                })
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
  reduce_single(.data, ...,
                rfun = function(x, ...) {
                  disjoin(x, ..., ignore.strand = FALSE)
                })
}

#' @method disjoin_ranges_directed GroupedGenomicRanges
#' @export
disjoin_ranges_directed.GroupedGenomicRanges <- function(.data, ...) {
  reduce_by_grp(.data, ...,
                rfun = function(x, ...) {
                  disjoin(x, ..., ignore.strand = FALSE)
                })
}
