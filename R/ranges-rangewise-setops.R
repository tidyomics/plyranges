# ranges-rangewise-setops.R

#' Vector-wise Range set-operations
#'
#' @param x,y Two Ranges objects to compare.
#'
#' @details These are usual set-operations that act on the sets of the
#' ranges represented in x and y. By default these operations will ignore
#' any strand information. The directed versions of these functions will
#' take into account strand for GRanges objects.
#'
#' @return A Ranges object
#'
#' @examples
#' gr1 <- data.frame(seqnames = "chr1",
#'                   start = c(2,9),
#'                   end = c(7,9),
#'                   strand = c("+", "-")) %>%
#'                as_granges()
#' gr2 <- data.frame(seqnames = "chr1", start = 5, width = 5, strand = "-") %>%
#'          as_granges()
#'
#' union_ranges(gr1, gr2)
#' union_ranges_directed(gr1, gr2)
#'
#' intersect_ranges(gr1, gr2)
#' intersect_ranges_directed(gr1, gr2)
#'
#' setdiff_ranges(gr1, gr2)
#' setdiff_ranges_directed(gr1, gr2)
#' # taking the complement of a ranges requires annotation information
#' gr1 <- set_genome_info(gr1, seqlengths = 100)
#' complement_ranges(gr1)
#'
#' @rdname ranges-setops
#' @export
intersect_ranges <- function(x,y) { UseMethod("intersect_ranges") }

#' @importFrom IRanges intersect
#' @export
intersect_ranges.IntegerRanges <- function(x,y) {
  intersect(x,y)
}

#' @importFrom GenomicRanges intersect
#' @export
intersect_ranges.GenomicRanges <- function(x,y) {
  intersect(x,y,ignore.strand = TRUE)
}

#' @export
#' @rdname ranges-setops
intersect_ranges_directed <- function(x,y) { UseMethod("intersect_ranges_directed") }

#' @export
#' @importFrom GenomicRanges intersect
intersect_ranges_directed.GenomicRanges <- function(x,y) {
  intersect(x,y, ignore.strand = FALSE)
}

#' @export
#' @rdname ranges-setops
union_ranges <- function(x,y) { UseMethod("union_ranges") }

#' @export
#' @importFrom IRanges union
union_ranges.IntegerRanges <- function(x,y) {
  union(x,y)
}

#' @export
#' @importFrom GenomicRanges union
union_ranges.GenomicRanges <- function(x,y) {
  union(x,y,ignore.strand = TRUE)
}

#' @export
#' @rdname ranges-setops
union_ranges_directed <- function(x,y) { UseMethod("union_ranges_directed") }

#' @export
#' @importFrom GenomicRanges union
union_ranges_directed.GenomicRanges <- function(x,y) {
  union(x,y, ignore.strand = FALSE)
}

#' @export
#' @rdname ranges-setops
setdiff_ranges <- function(x,y) { UseMethod("setdiff_ranges") }

#' @export
#' @importFrom IRanges setdiff
setdiff_ranges.IntegerRanges <- function(x,y) {
  setdiff(x,y)
}

#' @export
#' @importFrom GenomicRanges setdiff
setdiff_ranges.GenomicRanges <- function(x,y) {
  setdiff(x,y, ignore.strand = TRUE)
}

#' @export
#' @rdname ranges-setops
setdiff_ranges_directed <- function(x,y) { UseMethod("setdiff_ranges_directed") }

#' @export
#' @importFrom GenomicRanges setdiff
setdiff_ranges_directed.GenomicRanges <- function(x,y) {
  setdiff(x,y, ignore.strand = FALSE)
}

#' @rdname ranges-setops
#' @export
complement_ranges  <- function(x) { UseMethod("complement_ranges") }

#' @export
complement_ranges.IntegerRanges <- function(x) {
  setdiff(IRanges(start = min(start(x)), end = max(end(x))), x)
}

#' @export
complement_ranges.GenomicRanges <- function(x) {
  y <- try(get_genome_info(x), silent = TRUE)
  if (is(y, "try-error")) {
    stop("complement_ranges requires sequence lengths information", call. = FALSE)
  }
  setdiff_ranges(y, x)
}

#' @rdname ranges-setops
#' @export
complement_ranges_directed <- function(x) {UseMethod("complement_ranges_directed")}

#' @export
complement_ranges_directed.GenomicRanges <- function(x) {
  y <- try(get_genome_info(x), silent = TRUE)
  if (is(y, "try-error")) {
    stop("complement_ranges_directed requires sequence lengths information", call. = FALSE)
  }
  setdiff_ranges_directed(y, x)
}