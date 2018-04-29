#' Functional setters for Ranges objects
#' @details These methods are used internally in `mutate()` to modify 
#' core columns in Ranges objects.
#' @param x a Ranges object
#' @param width integer amount to modify width by
#' @param start integer amount to modify start by
#' @param end integer amount to modify end by
#' @param seqnames update seqnames column
#' @param strand update strand column
#' 
#' @rdname ranges-setters
#' @importFrom IRanges resize
#' @export
set_width <- function(x, width) UseMethod("set_width")

set_width.Ranges <- function(x, width = 0L) {
  width(x) <- width
  x
}

set_width.AnchoredIntegerRanges <- function(x, width = 0L) {
  anchor <- anchor(x)
  rng <- x@delegate
  switch(anchor,
         start = resize(rng, width, fix = "start"),
         end = resize(rng, width, fix = "end"),
         center = resize(rng, width, fix = "center")
  )
}


set_width.AnchoredGenomicRanges <- function(x, width = 0L) {
  anchor <- anchor(x)
  rng <- x@delegate
  switch(anchor,
         start = resize(rng, width, fix = "start", ignore.strand = TRUE),
         end = resize(rng, width, fix = "end", ignore.strand = TRUE),
         center = resize(rng, width, fix = "center", ignore.strand = TRUE),
         "3p" = resize_by_strand(rng, width, "3p"),
         "5p" = resize_by_strand(rng, width, "5p")
  )
}

# anchor by strand version of resize
resize_by_strand <- function(x, width, anchor) {
  # anchor_3p == fix by start for negative strand
  if (anchor == "3p") {
    return(resize(x, width, fix = "end", ignore.strand = FALSE))
  } else {
    return(resize(x, width, fix = "start", ignore.strand = FALSE))
  }
}

#' @rdname ranges-setters
#' @export
set_start <- function(x, start = 0L) { UseMethod("set_start") }

set_start.Ranges <- function(x, start = 0L) {
  start(x) <- start
  x
}

#' @rdname ranges-setters
#' @export
set_end <- function(x, end = 0L) { UseMethod("set_end") }

set_end.Ranges <- function(x, end = 0L) {
  end(x) <- end
  x
}

#' @rdname ranges-setters
#' @export
set_seqnames <- function(x, seqnames) { UseMethod("set_seqnames") }

set_seqnames.GenomicRanges <- function(x, seqnames) {
  seqnames(x) <- seqnames
  x
}

#' @rdname ranges-setters
#' @export
set_strand <- function(x, strand) { UseMethod("set_strand") }

set_strand.GenomicRanges <- function(x, strand = "*") {
  strand(x) <- strand
  x
}