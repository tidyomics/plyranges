#' Internal methods for `mutate`
#' functional setters for Ranges objects
set_start <- function(x, start = 0L) { UseMethod("set_start") }

set_start.Ranges <- function(x, start = 0L) {
  start(x) <- start
  x
}

set_end <- function(x, end = 0L) { UseMethod("set_start") }

set_end.Ranges <- function(x, end = 0L) {
  end(x) <- end
  x
}

set_seqnames <- function(x, seqnames) { UseMethod("set_seqnames") }

set_seqnames.GenomicRanges <- function(x, seqnames) {
  seqnames(x) <- seqnames
  x
}

set_strand <- function(x, seqnames) { UseMethod("set_seqnames") }

set_strand.GenomicRanges <- function(x, strand = "*") {
  strand(x) <- strand
  x
}