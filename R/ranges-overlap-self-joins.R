#' Find overlaps within a Ranges object
#'
#' @param x A Ranges object
#' @param maxgap,minoverlap The maximimum gap between intervals as an integer
#' greater than or equal to zero. The minimum amount of overlap between intervals
#' as an integer greater than zero, accounting for the maximum gap.
#'
#' @details Self overlaps find any overlaps (or overlaps within or overlaps
#' directed) between a ranges object and itself.
#'
#' @return a Ranges object
#'
#' @examples
#' query <- data.frame(start = c(5,10, 15,20), width = 5, gc = runif(4)) %>%
#'              as_iranges()
#'
#' join_overlap_self(query)

#'
#' # -- GRanges objects, strand is ignored by default
#' query  <- data.frame(seqnames = "chr1",
#'                start = c(11,101),
#'                end = c(21, 200),
#'                name = c("a1", "a2"),
#'                strand = c("+", "-"),
#'                score = c(1,2)) %>%
#'            as_granges()
#'
#' # ignores strandedness
#' join_overlap_self(query)
#' join_overlap_self_within(query)
#' # adding directed prefix includes strand
#' join_overlap_self_directed(query)
#'
#'
#' @rdname ranges-overlaps-self
#' @seealso [find_overlaps()][join_overlap_inner()]
#' @export
join_overlap_self <- function(x, maxgap, minoverlap) {
  UseMethod("join_overlap_self")
}

#' @export
join_overlap_self.GenomicRanges <- function(x, maxgap = -1L, minoverlap = 0L) {
  find_overlaps(x,x, maxgap, minoverlap, suffix = c("", ".overlap"))
}

#' @export
join_overlap_self.IntegerRanges <- join_overlap_self.GenomicRanges


#' @rdname ranges-overlaps-self
#' @export
join_overlap_self_within <- function(x, maxgap, minoverlap) {
  UseMethod("join_overlap_self_within")
}

#' @export
join_overlap_self_within.GenomicRanges <- function(x, maxgap = -1L, minoverlap = 0L) {
  find_overlaps_within(x,x, maxgap, minoverlap, suffix = c("", ".overlap"))
}

#' @export
join_overlap_self_within.IntegerRanges <- join_overlap_self_within.GenomicRanges

#' @rdname ranges-overlaps-self
#' @export
join_overlap_self_directed <- function(x, maxgap, minoverlap) {
  UseMethod("join_overlap_self_directed")
}

#' @export
join_overlap_self_directed.GenomicRanges <- function(x, maxgap = -1L, minoverlap = 0L) {
  find_overlaps_directed(x,x, maxgap, minoverlap, suffix = c("", ".overlap"))
}