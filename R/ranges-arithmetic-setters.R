#' Internal methods for `mutate`
#' functional setters for Ranges objects
#' @export
set_start <- function(x, start = 0L) { UseMethod("set_start") }

set_start.Ranges <- function(x, start = 0L) {
  start(x) <- start
  x
}

#' @export
set_end <- function(x, end = 0L) { UseMethod("set_end") }

set_end.Ranges <- function(x, end = 0L) {
  end(x) <- end
  x
}

#' @export
set_seqnames <- function(x, seqnames) { UseMethod("set_seqnames") }

set_seqnames.GenomicRanges <- function(x, seqnames) {
  seqnames(x) <- seqnames
  x
}

#' @export
set_strand <- function(x, strand) { UseMethod("set_strand") }

set_strand.GenomicRanges <- function(x, strand = "*") {
  strand(x) <- strand
  x
}