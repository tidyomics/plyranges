#' @rdname ranges-overlaps
#' @importFrom rlang syms
#' @export
group_by_overlaps <- function(x, y, maxgap, minoverlap) { UseMethod("group_by_overlaps") }

#' @rdname ranges-overlaps
#' @export
group_by_overlaps.IntegerRanges <- function(x, y, maxgap = -1L, minoverlap = 0L) {
  hits <- make_hits(x, y, findOverlaps, maxgap = maxgap, minoverlap = minoverlap)
  left <- expand_by_hits(x, y, c(".query", ".subject"), hits)
  mcols(left)$query <- queryHits(hits)
  new_grouped_ir(left, !!!rlang::syms("query"))
}

#' @rdname ranges-overlaps
#' @export
group_by_overlaps.GenomicRanges <- function(x, y, maxgap = -1L, minoverlap = 0L) {
  hits <- make_hits(x, y, findOverlaps, 
                    maxgap = maxgap, 
                    minoverlap = minoverlap,
                    ignore.strand = TRUE)
  left <- expand_by_hits(x, y, c(".query", ".subject"), hits)
  mcols(left)$query <- queryHits(hits)
  new_grouped_gr(left, !!!rlang::syms("query"))
}

# TODO add in more variants here?