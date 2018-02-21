# ranges-reduce
reduce_rng <- function(.data, reduced, dots) {

  revmap <- mcols(reduced)$revmap

  ranges_list <- IRanges::relist(.data[unlist(revmap)], revmap)

  reduced_summary <- as(lapply(ranges_list, summarize_rng, dots), "List")

  mcols(reduced) <- do.call(rbind, lapply(reduced_summary, as, "DataFrame"))
  return(reduced)
}

#' Reduce then aggregate a Ranges object
#'
#' @param .data a Ranges object to reduce
#' @param ... Name-value pairs of summary functions.
#'
#' @return a Ranges object with the
#' @rdname ranges-reduce
#' @importFrom IRanges reduce
#' @importFrom utils relist
#' @examples
#' df <- data.frame(start = 1:10, width = 5,  seqnames = "seq1",
#' strand = sample(c("+", "-", "*"), 10, replace = TRUE), gc = runif(10))
#' rng <- as_granges(df)
#' rng %>% reduce_ranges()
#' rng %>% reduce_ranges(gc = mean(gc))
#' rng %>% reduce_ranges_directed(gc = mean(gc))
#'
#' x <- data.frame(start = c(11:13, 2, 7:6), width=3, id=letters[1:6],
#' score=1:6)
#' x <- as_iranges(x)
#' x %>% reduce_ranges()
#' x %>% reduce_ranges(score = sum(score))
#' @export
reduce_ranges <- function(.data, ...) { UseMethod("reduce_ranges") }

#' @method reduce_ranges Ranges
#' @export
reduce_ranges.Ranges <- function(.data, ...) {
  dots <- quos(...)
  if (length(dots) == 0L) {
    return(reduce(.data))
  }

  reduced <- reduce(.data, with.revmap = TRUE)

  reduce_rng(.data, reduced, dots)

}


#' @method reduce_ranges GRangesGrouped
#' @export
reduce_ranges.GRangesGrouped <- function(.data, ...) {
  dots <- quos(...)
  split_ranges <- split_groups(.data, populate_mcols = TRUE, drop = TRUE)
  if (length(dots) == 0L) {
    gr_r <- IRanges::stack(reduce(split_ranges, ignore.strand = TRUE))
    mcols(gr_r) <- mcols(gr_r)[, group_vars(.data), drop = FALSE]
    return(gr_r)
  }

  reduced <- reduce(split_ranges, with.revmap = TRUE, ignore.strand = TRUE)

  gr_r <- S4Vectors::List(Map(function(i)
    reduce_rng(split_ranges[[i]], reduced[[i]], dots),
    seq_along(split_ranges)))
  mcols(gr_r) <- mcols(split_ranges)
  gr_r <- IRanges::stack(gr_r)
  mcols(gr_r) <- mcols(gr_r)[, -1, drop = FALSE]
  mcols(gr_r) <- mcols(gr_r)[, rev(seq_along(mcols(gr_r)))]
  return(gr_r)
}

#' @method reduce_ranges GenomicRanges
#' @export
reduce_ranges.GenomicRanges <- function(.data, ...) {

  dots <- quos(...)
  if (length(dots) == 0L) {
    return(reduce(.data,ignore.strand = TRUE))
  }

  reduced <- reduce(.data, with.revmap = TRUE, ignore.strand = TRUE)
  reduce_rng(.data, reduced, dots)
}

#' @rdname ranges-reduce
#' @export
reduce_ranges_directed <- function(.data, ...) {
  UseMethod("reduce_ranges_directed")
}

#' @importFrom IRanges reduce
#' @method reduce_ranges_directed GenomicRanges
#' @export
reduce_ranges_directed.GenomicRanges <- function(.data, ...) {

  dots <- quos(...)
  if (length(dots) == 0L) {
    return(reduce(.data,ignore.strand = FALSE))
  }

  reduced <- reduce(.data, with.revmap = TRUE, ignore.strand = FALSE)
  reduce_rng(.data, reduced, dots)

}
