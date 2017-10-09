# ranges-joins
mcols_overlaps_inner_update <- function(left, right, suffix) {

  left_names <- names(mcols(left))
  right_names <- names(mcols(right))
  common_name <- intersect(left_names, right_names)
  lname_inx <- left_names %in% common_name
  rname_inx <- right_names %in% common_name
  if (any(lname_inx)) {
    names(mcols(left))[lname_inx] <- paste0(left_names[lname_inx], suffix[1])
  }

  if (any(rname_inx)) {
    names(mcols(right))[rname_inx] <- paste0(right_names[rname_inx], suffix[2])
  }

  additional_mcols <- DataFrame(width.x = width(left),
                                width.y = width(right))
  if (!is.null(mcols(left))) {
    additional_mcols <- cbind(additional_mcols, mcols(left))
  }

  if (!is.null(mcols(right))) {
    additional_mcols <- cbind(additional_mcols, mcols(right))
  }
  additional_mcols
}

#' Join by overlapping Ranges
#'
#' @param x,y Objects representing ranges
#' @param maxgap,minoverlap The maximimum gap between intervals as an integer
#' greater than or equal to zero. The minimum amount of overlap between intervals
#' as an integer greater than zero, accounting for the maximum gap.
#' @param suffix Character to vectors to append to common columns in x and y.
#'
#' @details Inner overlaps finds intersecting intervals between x and y
#' while keeping information about the width from x and y.
#' Left overlaps is equivalent to find_overlaps, it keeps all ranges in x
#' that overlap ranges in y.
#'
#' @importFrom S4Vectors first second DataFrame
#' @importFrom GenomicRanges findOverlapPairs
#' @rdname overlap-joins
#' @export
join_overlap_inner <- function(x, y, maxgap, minoverlap, suffix = c(".x", ".y")) UseMethod("join_overlap_inner")

join_overlap_inner.Ranges <- function(x, y, maxgap = 0L, minoverlap = 1L,
                                      suffix = c(".x", ".y")) {
  pairs <- findOverlapPairs(x, y,
                            type = "any",
                            maxgap = maxgap,
                            minoverlap = minoverlap)
  inner_rng <- pintersect(pairs)
  left_rng <- first(pairs)
  right_rng <- second(pairs)
  mcols(inner_rng) <- mcols_overlaps_inner_update(left_rng, right_rng, suffix)
  inner_rng

}


join_overlap_inner.GenomicRanges <- function(x, y, maxgap = 0L, minoverlap = 1L,
                                             suffix = c(".x", ".y")) {
  pairs <- findOverlapPairs(x, y,
                            type = "any",
                            maxgap = maxgap,
                            minoverlap = minoverlap,
                            ignore.strand = TRUE)
  inner_rng <- pintersect(pairs, ignore.strand = TRUE)

  left_rng <- first(pairs)
  right_rng <- second(pairs)
  mcols(inner_rng) <- mcols_overlaps_inner_update(left_rng, right_rng, suffix)
  inner_rng

}

#' @rdname overlap-joins
join_overlap_within_inner <- function(x, y, maxgap, minoverlap,
                                      suffix = c(".x", ".y")) {
  UseMethod("join_overlap_inner")
}

join_overlap_within_inner.Ranges <- function(x, y, maxgap, minoverlap,
                                             suffix = c(".x", ".y")) {
  pairs <- findOverlapPairs(x, y,
                            type = "within",
                            maxgap = maxgap,
                            minoverlap = minoverlap)
  inner_rng <- pintersect(pairs)
  left_rng <- first(pairs)
  right_rng <- second(pairs)
  mcols(inner_rng) <- mcols_overlaps_inner_update(left_rng, right_rng, suffix)
  inner_rng
}

join_overlap_within_inner.GenomicRanges <- function(x, y, maxgap, minoverlap,
                                             suffix = c(".x", ".y")) {
  pairs <- findOverlapPairs(x, y,
                            type = "within",
                            maxgap = maxgap,
                            minoverlap = minoverlap,
                            ignore.strand = TRUE)
  inner_rng <- pintersect(pairs)
  left_rng <- first(pairs)
  right_rng <- second(pairs)
  mcols(inner_rng) <- mcols_overlaps_inner_update(left_rng, right_rng, suffix)
  inner_rng
}

#' @rdname overlap-joins
#' @export
join_overlap_left <- function(x, y, maxgap, minoverlap, suffix = c(".x", ".y")) {
  UseMethod("join_overlap_left")
}

join_overlap_left.Ranges <- function(x, y, maxgap = 0L, minoverlap = 1L, suffix = c(".x", ".y")) {
  find_overlaps(x,y,maxgap, mingap,suffix)
}

join_overlap_left.GenomicRanges <- function(x, y, maxgap = 0L, minoverlap = 1L, suffix = c(".x", ".y")) {
  find_overlaps(x,y,maxgap, mingap,suffix)
}



