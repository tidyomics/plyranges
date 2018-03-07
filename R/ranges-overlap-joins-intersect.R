#' Join by overlapping Ranges
#'
#' @param x,y Objects representing ranges
#' @param maxgap,minoverlap The maximimum gap between intervals as an integer
#' greater than or equal to zero. The minimum amount of overlap between intervals
#' as an integer greater than zero, accounting for the maximum gap.
#' @param suffix Character to vectors to append to common columns in x and y
#' (default = \code{c(".x", ".y")}).
#'
#' @details The function \code{join_intersect_overlaps} finds
#' the genomic intervals that are the overlapping ranges between x and y and
#' returns a new ranges object with metadata columns from x and y.
#'
#' The function \code{join_inner_overlaps} is equivalent to \code{find_overlaps}
#' it returns all ranges in x that overlap ranges in y.
#'
#' The function \code{join_left_overlaps} performs a left outer join between x
#' and y. It returns all ranges in x that overlap or do not overlap ranges in y
#' plus metadata columns common to both. If there is no overlapping range
#' the metadata column will contain a missing value.
#'
#' The function \code{join_self_overlaps} find all overlaps between a ranges
#' object x and itself.
#'
#' All of these functions have two suffixes that modify their behavior.
#' The \code{within} suffix, returns only ranges in x that are completely
#' contained in y. The \code{directed} suffix takes into account the strandedness
#' of a GRanges object.
#'
#' @return a GRanges object
#'
#' @importFrom S4Vectors first second DataFrame
#' @importFrom IRanges findOverlapPairs
#' @examples
#' x <- as_iranges(data.frame(start = c(11, 101), end = c(21, 201)))
#' y <- as_iranges(data.frame(start = c(10, 20, 50, 100),
#'                            end = c(19, 21, 105, 202)))
#'
#' # intersect takes common interval
#' join_overlap_intersect(x,y)
#'
#' # within
#' join_overlap_intersect_within(x,y)
#'
#' @rdname overlap-joins
#' @export
#'
#'
join_overlap_intersect <- function(x, y, maxgap, minoverlap, suffix = c(".x", ".y")) {
  UseMethod("join_overlap_intersect")
}

#' @export
join_overlap_intersect.Ranges <- function(x, y, maxgap = -1L, minoverlap = 0L,
                                          suffix = c(".x", ".y")) {
  pairs <- findOverlapPairs(x, y,
                            type = "any",
                            maxgap = maxgap,
                            minoverlap = minoverlap)
  inner_rng <- pintersect(pairs)
  left_rng <- first(pairs)
  right_rng <- second(pairs)
  mcols(inner_rng) <- mcols_overlaps_update(left_rng, right_rng, suffix)
  inner_rng

}

#' @export
join_overlap_intersect.GenomicRanges <- function(x, y, maxgap = -1L, minoverlap = 0L,
                                                 suffix = c(".x", ".y")) {
  pairs <- findOverlapPairs(x, y,
                            type = "any",
                            maxgap = maxgap,
                            minoverlap = minoverlap,
                            ignore.strand = TRUE)

  inner_rng <- pintersect(pairs, ignore.strand = TRUE)

  left_rng <- first(pairs)
  right_rng <- second(pairs)
  mcols(inner_rng) <- mcols_overlaps_update(left_rng, right_rng, suffix)
  inner_rng

}
#' @rdname overlap-joins
#' @export
join_overlap_intersect_within <- function(x, y, maxgap, minoverlap,
                                      suffix = c(".x", ".y")) {
  UseMethod("join_overlap_intersect_within")
}

#' @export
join_overlap_intersect_within.Ranges <- function(x, y, maxgap = -1L, minoverlap = 0L,
                                             suffix = c(".x", ".y")) {
  pairs <- findOverlapPairs(x, y,
                            type = "within",
                            maxgap = maxgap,
                            minoverlap = minoverlap)
  inner_rng <- pintersect(pairs)
  left_rng <- first(pairs)
  right_rng <- second(pairs)
  mcols(inner_rng) <- mcols_overlaps_update(left_rng, right_rng, suffix)
  inner_rng
}

#' @export
join_overlap_intersect_within.GenomicRanges <- function(x, y, maxgap = -1L, minoverlap = 0L,
                                                    suffix = c(".x", ".y")) {
  pairs <- findOverlapPairs(x, y,
                            type = "within",
                            maxgap = maxgap,
                            minoverlap = minoverlap,
                            ignore.strand = TRUE)
  inner_rng <- pintersect(pairs)
  left_rng <- first(pairs)
  right_rng <- second(pairs)
  mcols(inner_rng) <- mcols_overlaps_update(left_rng, right_rng, suffix)
  inner_rng
}


#' @rdname overlap-joins
#' @export
join_overlap_intersect_directed <- function(x, y, maxgap, minoverlap,
                                            suffix = c(".x", ".y")) {
  UseMethod("join_overlap_intersect_directed")
}

#' @export
join_overlap_intersect_directed.GenomicRanges <- function(x, y,
                                                          maxgap = -1L,
                                                          minoverlap = 0L,
                                                          suffix = c(".x", ".y")) {
  pairs <- findOverlapPairs(x, y,
                            type = "any",
                            maxgap = maxgap,
                            minoverlap = minoverlap,
                            ignore.strand = FALSE)

  inner_rng <- pintersect(pairs, ignore.strand = FALSE)

  left_rng <- first(pairs)
  right_rng <- second(pairs)
  mcols(inner_rng) <- mcols_overlaps_update(left_rng, right_rng, suffix)
  inner_rng
}
