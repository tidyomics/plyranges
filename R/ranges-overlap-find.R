# ranges-overlaps.R
# workhorse funciton for copying metadata columns in y
mcols_overlaps_update <- function(x, y, hits, suffix) {

  left <- x[queryHits(hits), ]
  right <- y[subjectHits(hits), ]
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

  additional_cols <- DataFrame(start = start(right),
                               end = end(right),
                               width = width(right))
  if (is(right, "GRanges")) additional_cols$strand <- strand(right)
  names(additional_cols) <- paste0(names(additional_cols), suffix[2])

  if (!is.null(mcols(left))) {
    additional_cols <- cbind(additional_cols, mcols(left))
  }

  if (!is.null(mcols(right))) {
    additional_cols <- cbind(additional_cols, mcols(right))
  }
  mcols(left) <- additional_cols

  left
}


#' Find overlap between two Ranges
#'
#' @rdname ranges-overlaps
#'
#' @param x,y Objects representing ranges
#' @param maxgap,minoverlap The maximimum gap between intervals as an integer
#' greater than or equal to negative one. The minimum amount of overlap between intervals
#' as an integer greater than zero, accounting for the maximum gap.
#' @param suffix A character vector length two used to identify metadata columns
#' coming from x and y
#'
#' @details \code{find_overlaps} will search for any overlaps between ranges
#' x and y and return a ranges object of the same length as x but with additional
#' metadata colums indicating the start, end, width of the overlap range in y
#' and any additional metadata columns in y. \code{find_overlaps_within} is
#' the same but will only search for overlaps within y. For GRanges strand is
#' ignored.
#'
#' @return A Ranges object with rows corresponding to the
#' ranges in x that overlap y. In the case of \code{group_by_overlaps}, returns
#' a GroupedRanges object, grouped by the number of overlaps
#' of ranges in x that overlap y (stored in a column called query).
#' @seealso \link[GenomicRanges]{setops-methods}, \link[IRanges]{findOverlaps-methods}
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @export
find_overlaps <- function(x, y, maxgap, minoverlap, suffix = c(".x", ".y")) {
  UseMethod("find_overlaps")
}

#' @rdname ranges-overlaps
#' @export
find_overlaps.Ranges <- function(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  hits <- findOverlaps(x,y, maxgap, minoverlap, type = "any", select = "all")
  mcols_overlaps_update(x,y, hits, suffix)
}

#' @rdname ranges-overlaps
#' @export
find_overlaps.GenomicRanges <- function(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  hits <- findOverlaps(x,y, maxgap, minoverlap,
                       type = "any", select = "all", ignore.strand = TRUE)
  mcols_overlaps_update(x,y,hits, suffix)
}

#' @rdname ranges-overlaps
#' @export
find_overlaps_within <- function(x, y, maxgap, minoverlap, suffix = c(".x", ".y")) {
  UseMethod("find_overlaps_within")
}

#' @rdname ranges-overlaps
#' @export
find_overlaps_within.Ranges <- function(x,y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  hits <- findOverlaps(x,y, maxgap, minoverlap, type = "within", select = "all")
  mcols_overlaps_update(x,y, hits, suffix = c(".x", ".y"))
}

#' @rdname ranges-overlaps
#' @export
find_overlaps_within.GenomicRanges <- function(x,y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  hits <- findOverlaps(x,y, maxgap, minoverlap, type = "within", select = "all")
  mcols_overlaps_update(x,y, hits, suffix = c(".x", ".y"))
}
