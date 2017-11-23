#' @rdname ranges-overlaps
#' @importFrom rlang syms
#' @export
group_by_overlaps <- function(x, y, maxgap, minoverlap) { UseMethod("group_by_overlaps") }

#' @rdname ranges-overlaps
#' @export
group_by_overlaps.Ranges <- function(x, y, maxgap = -1L, minoverlap = 0L) {

  hits <- findOverlaps(x, y, maxgap, minoverlap,
                       type = "any", select = "all")
  rng <- mcols_overlaps_update(x,y,hits, suffix = c(".query", ".subject"))
  mcols(rng)$query <- queryHits(hits)
  new("IRangesGrouped", rng, groups = syms("query"))
}

#' @rdname ranges-overlaps
#' @export
group_by_overlaps.GenomicRanges <- function(x, y, maxgap = -1L, minoverlap = 0L) {

  hits <- findOverlaps(x,y, maxgap, minoverlap,
                       type = "any", select = "all", ignore.strand = TRUE)
  rng <- mcols_overlaps_update(x,y,hits, suffix = c(".query", ".subject"))
  mcols(rng)$query <- queryHits(hits)
  new("GRangesGrouped", rng, groups = syms("query"))
}


