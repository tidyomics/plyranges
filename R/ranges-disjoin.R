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
#' rng <- as_granges(df)
#' rng %>% disjoin_ranges()
#' rng %>% disjoin_ranges(gc = mean(gc))
#' rng %>% disjoin_ranges_directed(gc = mean(gc))
#' @export
disjoin_ranges <- function(.data, ...) { UseMethod("disjoin_ranges") }

#' @method disjoin_ranges IntegerRanges
#' @export
disjoin_ranges.IntegerRanges <- function(.data, ...) {
  dots <- quos(...)
  if (length(dots) == 0L) {
    return(disjoin(.data))
  }

  disjoined <- disjoin(.data, with.revmap = TRUE)

  reduce_rng(.data, disjoined, dots)

}

#' @method disjoin_ranges GroupedIntegerRanges
#' @export
disjoin_ranges.GroupedIntegerRanges <- function(.data, ...) {
  dots <- quos(...)
  delegate <- .data@delegate
  inx <- .data@inx
  split_ranges <- extractList(delegate, inx)
  mcols(split_ranges) <- mcols(inx)
  
  if (length(dots) == 0L) {
    ir_d <- IRanges::stack(disjoin(split_ranges, ignore.strand = TRUE))
    mcols(ir_d) <- mcols(ir_d)[, group_vars(.data), drop = FALSE]
    return(ir_d)
  }

  disjoined <- disjoin(split_ranges, with.revmap = TRUE, ignore.strand = TRUE)

  ir_d <- S4Vectors::List(Map(function(i)
    reduce_rng(split_ranges[[i]], disjoined[[i]], dots),
    seq_along(split_ranges)))
  mcols(ir_d) <- mcols(split_ranges)
  ir_d <- IRanges::stack(ir_d)
  mcols(ir_d) <- mcols(ir_d)[, -1, drop = FALSE]
  mcols(ir_d) <- mcols(ir_d)[, rev(seq_along(mcols(ir_d)))]
  return(ir_d)
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

#' @method disjoin_ranges GroupedGenomicRanges
#' @export
disjoin_ranges.GroupedGenomicRanges <- function(.data, ...) {
  dots <- quos(...)
  delegate <- .data@delegate
  inx <- .data@inx
  split_ranges <- extractList(delegate, inx)
  mcols(split_ranges) <- mcols(inx)
  
  if (length(dots) == 0L) {
    gr_d <- IRanges::stack(disjoin(split_ranges, ignore.strand = TRUE))
    mcols(gr_d) <- mcols(gr_d)[, group_vars(.data), drop = FALSE]
    return(gr_d)
  }

  disjoined <- disjoin(split_ranges, with.revmap = TRUE, ignore.strand = TRUE)

  gr_d <- S4Vectors::List(Map(function(i)
    reduce_rng(split_ranges[[i]], disjoined[[i]], dots),
    seq_along(split_ranges)))
  mcols(gr_d) <- mcols(split_ranges)
  gr_d <- IRanges::stack(gr_d)
  mcols(gr_d) <- mcols(gr_d)[, -1, drop = FALSE]
  mcols(gr_d) <- mcols(gr_d)[, rev(seq_along(mcols(gr_d)))]
  return(gr_d)
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

#' @method disjoin_ranges_directed GroupedGenomicRanges
#' @export
disjoin_ranges_directed.GroupedGenomicRanges <- function(.data, ...) {
  dots <- quos(...)
  delegate <- .data@delegate
  inx <- .data@inx
  split_ranges <- extractList(delegate, inx)
  mcols(split_ranges) <- mcols(inx)
  
  if (length(dots) == 0L) {
    gr_d <- IRanges::stack(disjoin(split_ranges, ignore.strand = FALSE))
    mcols(gr_d) <- mcols(gr_d)[, group_vars(.data), drop = FALSE]
    return(gr_d)
  }

  disjoined <- disjoin(split_ranges, with.revmap = TRUE, ignore.strand = FALSE)

  gr_d <- S4Vectors::List(Map(function(i)
    reduce_rng(split_ranges[[i]], disjoined[[i]], dots),
    seq_along(split_ranges)))
  mcols(gr_d) <- mcols(split_ranges)
  gr_d <- IRanges::stack(gr_d)
  mcols(gr_d) <- mcols(gr_d)[, -1, drop = FALSE]
  mcols(gr_d) <- mcols(gr_d)[, rev(seq_along(mcols(gr_d)))]
  return(gr_d)
}
