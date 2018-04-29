#' Count the number of overlaps between two Ranges objects
#'
#' @param x,y Objects representing ranges
#' @param maxgap,minoverlap The maximimum gap between intervals as an integer
#' greater than or equal to zero. The minimum amount of overlap between intervals
#' as an integer greater than zero, accounting for the maximum gap.
#'
#' @return An integer vector of same length as x.
#'
#' @examples
#' query <- data.frame(start = c(5,10, 15,20), width = 5, gc = runif(4)) %>%
#'              as_iranges()
#' subject <- data.frame(start = 2:6, width = 3:7, label = letters[1:5]) %>%
#'              as_iranges()
#' query %>% mutate(n_olap = count_overlaps(., subject),
#'                  n_olap_within = count_overlaps_within(., subject))
#'
#' @importFrom IRanges countOverlaps
#' @rdname ranges-count-overlaps
#' @export
count_overlaps <- function(x, y, maxgap, minoverlap) {
  UseMethod("count_overlaps")
}

#' @rdname ranges-count-overlaps
#' @export
count_overlaps.IntegerRanges <- function(x,y, maxgap = -1L, minoverlap = 0L) {
  countOverlaps(x,y, maxgap, minoverlap, type = "any")
}

#' @rdname ranges-count-overlaps
#' @export
count_overlaps.GenomicRanges <- function(x,y, maxgap = -1L, minoverlap = 0L) {
  countOverlaps(x,y, maxgap, minoverlap, type = "any", ignore.strand = TRUE)
}

#' @rdname ranges-count-overlaps
#' @export
count_overlaps_within <- function(x, y, maxgap, minoverlap) {
  UseMethod("count_overlaps_within")
}

#' @rdname ranges-count-overlaps
#' @export
count_overlaps_within.IntegerRanges <- function(x,y, maxgap = 0L, minoverlap = 1L) {
  countOverlaps(x,y, maxgap, minoverlap, type = "within")
}

#' @rdname ranges-count-overlaps
#' @export
count_overlaps_within.GenomicRanges <- function(x,y, maxgap = 0L, minoverlap = 1L) {
  countOverlaps(x,y, maxgap, minoverlap, type = "within", ignore.strand = TRUE)
}

#' @rdname ranges-count-overlaps
#' @export
count_overlaps_directed <- function(x, y, maxgap, minoverlap) {
  UseMethod("count_overlaps_directed")
}

#' @rdname ranges-count-overlaps
#' @export
count_overlaps_directed.GenomicRanges <- function(x,y, maxgap = -1L, minoverlap = 0L) {
  countOverlaps(x,y, maxgap, minoverlap, type = "any", ignore.strand = FALSE)
}

#' @rdname ranges-count-overlaps
#' @export
count_overlaps_within_directed <- function(x, y, maxgap, minoverlap) {
  UseMethod("count_overlaps_directed")
}

#' @rdname ranges-count-overlaps
#' @export
count_overlaps_within_directed.GenomicRanges <- function(x,y, maxgap = -1L, minoverlap = 0L) {
  countOverlaps(x,y, maxgap, minoverlap, type = "within", ignore.strand = FALSE)
}