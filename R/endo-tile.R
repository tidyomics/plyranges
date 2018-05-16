partition_rng <- function(rglist) {
  n <- elementNROWS(rglist)
  partitions <- rep.int(seq_along(n), n)
  res <- rglist@unlistData
  mcols(res)[["partition"]] <- partitions
  res
}

#' Make windows over a Ranges object
#' 
#' @param x a Ranges object
#' @param width the maximum width of each window/tile
#' 
#' @importFrom IRanges tile
#' @export
tile_ranges <- function(x, width) { UseMethod("tile_ranges") }

#' @export
tile_ranges.Ranges <- function(x, width) {
  rng <- tile(x, width = width)
  partition_rng(rng)
}

#' Make sliding windows over a Ranges object
#' 
#' @param x a Ranges object
#' @param width the maximum width of each window/tile
#' @param step the distance between start position of each sliding window 
#' 
#' @importFrom IRanges slidingWindows
#' @export
slide_ranges <- function(x, width) { UseMethod("tile_ranges") }

#' @export
slide_ranges.Ranges <- function(x, width, step) {
  rng <- slidingWindows(x, width = width, step = step)
  partition_rng(rng)
}