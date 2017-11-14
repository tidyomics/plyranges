
#' Group by overlaps
#'
#' @param x,y Objects representing ranges
#' @param maxgap,minoverlap The maximimum gap between intervals as an integer
#' greater than or equal to zero. The minimum amount of overlap between intervals
#' as an integer greater than zero, accounting for the maximum gap.
#'
#' @rdname ranges-overlaps.Rd
#' @export
group_by_overlaps <- function(x, y, maxgap, minoverlap) { UseMethod("group_by_overlaps") }

#' @rdname ranges-overlaps.Rd
#' @export
group_by_overlaps.Ranges <- function(x, y, maxgap = -1L, minoverlap = 0L) {

  hits <- findOverlaps(x, y, maxgap, minoverlap,
                       type = "any", select = "all")
  rng <- mcols_overlaps_update(x,y,hits, suffix = c(".query", ".subject"))
  mcols(rng)$query <- queryHits(hits)
  new("IRangesGrouped", rng, groups = syms("query"))
}

#' @rdname ranges-overlaps.Rd
#' @export
group_by_overlaps.GenomicRanges <- function(x, y, maxgap = -1L, minoverlap = 0L) {

  hits <- findOverlaps(x,y, maxgap, minoverlap,
                       type = "any", select = "all", ignore.strand = TRUE)
  rng <- mcols_overlaps_update(x,y,hits, suffix = c(".query", ".subject"))
  mcols(rng)$query <- queryHits(hits)
  new("IRangesGrouped", rng, groups = syms("query"))
}


