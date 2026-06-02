partition_rng <- function(rglist) {
  n <- elementNROWS(rglist)
  partitions <- rep.int(seq_along(n), n)
  res <- rglist@unlistData
  mcols(res)[["partition"]] <- partitions
  res
}

rearrange_partition_by_strand <- function(res) {
  if (is.null(mcols(res)[["partition"]])) {
    stop("'partition' metadata column is missing from ranges object")
  }
  n <- length(res)
  strands <- as.character(strand(res))
  dir <- ifelse(strands == "-", -1L, 1L)
  res[order(mcols(res)[["partition"]], dir * seq_len(n))]
}


#' @rdname ranges-tile
#' @importFrom IRanges tile
#' @export
tile_ranges <- function(x, width) { UseMethod("tile_ranges") }

#' @export
tile_ranges.Ranges <- function(x, width) {
  rng <- tile(x, width = width)
  partition_rng(rng)
}

#' @rdname ranges-tile
#' @export
tile_ranges_directed <- function(x, width) { UseMethod("tile_ranges_directed") }

#' @export
tile_ranges_directed.GenomicRanges <- function(x, width) {
  rearrange_partition_by_strand(tile_ranges(x, width))
}

#' Slide or tile over a Ranges object
#'
#' @param x a Ranges object
#' @param width the maximum width of each window/tile (integer vector of length 1)
#' @param step the distance between start position of each sliding window (integer vector of length 1)
#'
#' @details The `tile_ranges()` function paritions a Ranges object `x` by the given the
#' `width` over all ranges in `x`, truncated by the sequence end.
#' The `slide_ranges()` function makes sliding windows within each range of `x`
#' of size `width` and sliding by `step`.
#' Both `slide_ranges()` and `tile_ranges()` return a new Ranges object
#' with a metadata column called "partition" which contains the index of the
#' input range `x` that a parition belongs to.
#' The `_directed` variants order tiles/windows on negative strand ranges from
#' right to left (high to low coordinates), matching the direction of
#' transcription. They only accept `GRanges` inputs.
#'
#' @return a Ranges object
#' @importFrom IRanges slidingWindows
#' @examples
#'
#' gr <- data.frame(seqnames = c("chr1", rep("chr2", 3), rep("chr1", 2), rep("chr3", 4)),
#'                  start = 1:10,
#'                  end = 11,
#'                  strand = c("-", rep("+", 2), rep("*", 2), rep("+", 3), rep("-", 2))) %>%
#'       as_granges() %>%
#'       set_genome_info(seqlengths = c(11,12,13))
#'
#' # partition ranges into subranges of width 2, odd width ranges
#' # will have one subrange of width 1
#' tile_ranges(gr, width = 2)
#'
#' # directed: negative strand tiles ordered right to left
#' tile_ranges_directed(gr, width = 2)
#'
#' # make sliding windows of width 3, moving window with step size of 2
#' slide_ranges(gr, width = 3, step = 2)
#'
#' # directed: negative strand windows ordered right to left
#' slide_ranges_directed(gr, width = 3, step = 2)
#'
#' @seealso \code{GenomicRanges::\link[GenomicRanges:tile-methods]{tile()}}
#' @rdname ranges-tile
#' @export
slide_ranges <- function(x, width, step) { UseMethod("slide_ranges") }

#' @export
slide_ranges.Ranges <- function(x, width, step) {
  rng <- slidingWindows(x, width = width, step = step)
  partition_rng(rng)
}

#' @rdname ranges-tile
#' @export
slide_ranges_directed <- function(x, width, step) { UseMethod("slide_ranges_directed") }

#' @export
slide_ranges_directed.GenomicRanges <- function(x, width, step) {
  rearrange_partition_by_strand(slide_ranges(x, width, step))
}