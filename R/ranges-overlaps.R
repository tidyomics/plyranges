# ranges-overlaps.R

# workhorse funciton for copying metadata columns in y
# probably should be deleted once the left join is implemented.
mcols_overlaps_update <- function(x, y, hits) {
  hits_y <- subjectHits(hits)
  hits_x <- queryHits(hits)
  expand_y <- y[hits_y]
  # update start, end, width columns
  mcols(x)[["start.y"]] <- NA_integer_
  mcols(x)[["start.y"]][hits_x] <- start(expand_y)
  mcols(x)[["end.y"]] <- NA_integer_
  mcols(x)[["end.y"]][hits_x] <- end(expand_y)
  mcols(x)[["width.y"]] <- NA_integer_
  mcols(x)[["width.y"]][hits_x] <- width(expand_y)
  # copy metadata columns from y to x
  cols_to_add <- names(mcols(expand_y))
  for (col in cols_to_add) {
    cls <- class(mcols(expand_y)[[col]])
    if (col %in% names(mcols(x))) {
      col_x <- paste0(col, ".y")
    } else {
      col_x <- col
    }
    # this won't work if the column is not an S3 class or an atomic vector
    mcols(x)[[col_x]] <- vector(cls, length(x))
    mcols(x)[[col_x]][hits_x] <- mcols(expand_y)[[col]]
    mcols(x)[[col_x]][-hits_x] <- NA
  }
  return(x)
}


#' Find overlap between two Ranges
#'
#' @rdname ranges-overlaps.Rd
#'
#' @param x,y Objects representing ranges
#' @param maxgap,minoverlap The maximimum gap between intervals as an integer
#' greater than or equal to zero. The minimum amount of overlap between intervals
#' as an integer greater than zero, accounting for the maximum gap.
#'
#' @details \code{find_overlaps} will search for any overlaps between ranges
#' x and y and return a ranges object of the same length as x but with additional
#' metadata colums indicating the start, end, width of the overlap range in y
#' and any additional metadata columns in y. \code{find_overlaps_within} is
#' the same but will only search for overlaps within y. For GRanges strand is
#' ignored.
#'
#' @seealso \link[GenomicRanges]{setops-methods}, \link[IRanges]{findOverlaps-methods}
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @export
find_overlaps <- function(x, y, maxgap, minoverlap) {
  UseMethod("find_overlaps")
}

#' @rdname ranges-overlaps.Rd
#' @export
find_overlaps.Ranges <- function(x, y, maxgap = 0L, minoverlap = 1L) {
  hits <- findOverlaps(x,y, maxgap, minoverlap, type = "any", select = "all")
  mcols_overlaps_update(x,y, hits)
}

#' @rdname ranges-overlaps.Rd
#' @export
find_overlaps.GenomicRanges <- function(x, y, maxgap = 0L, minoverlap = 1L) {
  hits <- findOverlaps(x,y, maxgap, minoverlap,
                       type = "any", select = "all", ignore.strand = TRUE)
  mcols_overlaps_update(x,y,hits)
}

#' @rdname ranges-overlaps.Rd
#' @export
find_overlaps_within <- function(x, y, maxgap, minoverlap) {
  UseMethod("find_overlaps_within")
}

#' @rdname ranges-overlaps.Rd
#' @export
find_overlaps_within.Ranges <- function(x,y, maxgap = 0L, minoverlap = 1L) {
  hits <- findOverlaps(x,y, maxgap, minoverlap, type = "within", select = "all")
  mcols_overlaps_update(x,y, hits)
}

#' @rdname ranges-overlaps.Rd
#' @export
find_overlaps_within.GenomicRanges <- function(x,y, maxgap = 0L, minoverlap = 1L) {
  hits <- findOverlaps(x,y, maxgap, minoverlap, type = "within", select = "all")
  mcols_overlaps_update(x,y, hits)
}
