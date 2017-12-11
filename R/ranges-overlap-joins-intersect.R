#' Join by overlapping Ranges
#'
#' @param x,y Objects representing ranges
#' @param maxgap,minoverlap The maximimum gap between intervals as an integer
#' greater than or equal to zero. The minimum amount of overlap between intervals
#' as an integer greater than zero, accounting for the maximum gap.
#' @param suffix Character to vectors to append to common columns in x and y.
#'
#' @details The function \code{join_intersect_overlaps} finds
#' the genomic intervals that are the overlapping ranges between x and y and
#' returns a new ranges object with metadata columns from x and y.
#'
#' The function \code{join_inner_overlaps} is equivalent to \code{find_overlaps}
#' , it returns all ranges in x that overlap ranges in y.
#'
#' The function \code{join_left_overlaps} performs a left outer join between x
#' and y. It returns all ranges in x that overlap or do not overlap ranges in y
#' plus metadata columns common to both. If there is no overlapping range
#' the metadata column will contain a missing value.
#'
#' The function \code{join_self_overlaps} find all overlaps between a ranges
#' object x and itself.
#'
#' @importFrom S4Vectors first second DataFrame
#' @importFrom IRanges findOverlapPairs
#' @rdname overlap-joins
#' @export
join_overlap_intersect <- function(x, y, maxgap, minoverlap, suffix = c(".x", ".y")) {
  UseMethod("join_overlap_intersect")
}

#' @rdname overlap-joins
#' @export
join_overlap_intersect.Ranges <- function(x, y, maxgap = -1L, minoverlap = 1L,
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

#' @rdname overlap-joins
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

#' @rdname overlap-joins
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
