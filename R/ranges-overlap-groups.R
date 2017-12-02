#' @rdname ranges-overlaps
#' @importFrom rlang syms
#' @export
group_by_overlaps <- function(x, y, maxgap, minoverlap) { UseMethod("group_by_overlaps") }

#' @rdname ranges-overlaps
#' @export
group_by_overlaps.Ranges <- function(x, y, maxgap = -1L, minoverlap = 0L) {

  hits <- findOverlaps(x, y, maxgap, minoverlap,
                       type = "any", select = "all")

  left <- x[queryHits(hits), ]
  right <- y[subjectHits(hits), ]
  mcols(left) <- mcols_overlaps_update(left, right,
                                       suffix = c(".query", ".subject"),
                                       copy_left = FALSE)
  mcols(left)$query <- queryHits(hits)
  new("IRangesGrouped", left, groups = syms("query"))
}

#' @rdname ranges-overlaps
#' @export
group_by_overlaps.GenomicRanges <- function(x, y, maxgap = -1L, minoverlap = 0L) {

  hits <- findOverlaps(x,y, maxgap, minoverlap,
                       type = "any", select = "all", ignore.strand = TRUE)
  left <- x[queryHits(hits), ]
  right <- y[subjectHits(hits), ]
  mcols(left) <- mcols_overlaps_update(left, right,
                                       suffix = c(".query", ".subject"),
                                       copy_left = FALSE)
  mcols(left)$query <- queryHits(hits)
  new("GRangesGrouped", left,  groups = syms("query"))
}


