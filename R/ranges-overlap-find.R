# workhorse function for copying metadata columns in y
mcols_overlaps_update <- function(left, right, suffix, return_data_frame = FALSE) {

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

  if (!is.null(mcols(left))) {
    additional_mcols <- mcols(left)
  } else {
    additional_mcols <- NULL
  }

  if (!is.null(mcols(right))) {
    if (is.null(additional_mcols)) {
      additional_mcols <- mcols(right)
    } else {
      additional_mcols <- cbind(additional_mcols, mcols(right))
    }
  }

  if (return_data_frame) {
    if (is(left, "GenomicRanges")) {
      ranges_df <- DataFrame(granges.x = granges(left),
                             granges.y = granges(right))
    } else {
      ranges_df <- DataFrame(ranges.x = ranges(left),
                             ranges.y = ranges(right))
    }
    names(ranges_df) <- paste0(gsub("\\..*", "" , names(ranges_df)), suffix)
    if (!is.null(additional_mcols)) {
      additional_mcols <- cbind(ranges_df, additional_mcols)
    } else {
      return(ranges_df)
    }
  }

  additional_mcols
}

# backend function for olaps and related methods
# takes the pair of ranges, a suffix, optional arguments
# and f_in an anyonymous function that generates a Hits object
make_hits <- function(x, y, f_in, ...) {
  f_in(x, y, ...)
}

# expand x,y, based on Hits object
expand_by_hits <- function(x, y, suffix, hits, return_data_frame = FALSE) {
  if (is(hits, "Hits")) {
    left <- x[queryHits(hits), ]
    right <- y[subjectHits(hits), ]
  } else {
    # for vector case
    no_hits_id <- !is.na(hits)
    left <- x[no_hits_id, ]
    right <- y[hits[no_hits_id], ]
  }

  if (return_data_frame) {
    return(mcols_overlaps_update(left, right, suffix, return_data_frame))
  }
  mcols(left) <-  mcols_overlaps_update(left, right, suffix)
  left
}

.find_overlaps <- function(x,y, suffix, f_in, ...) {
  hits <- make_hits(x, y, f_in, ...)
  expand_by_hits(x,y, suffix, hits)
}

#' Find overlap between two Ranges
#'
#' @rdname ranges-overlaps
#'
#' @param x,y Objects representing ranges
#' @param maxgap,minoverlap The maximimum gap between intervals as an integer
#' greater than or equal to negative one. The minimum amount of overlap between intervals
#' as an integer greater than zero, accounting for the maximum gap.
#' @param suffix A character vector of length two used to identify metadata columns
#' coming from x and y.
#'
#' @details `find_overlaps()` will search for any overlaps between ranges
#' x and y and return a Ranges object of length equal to the number of times x
#' overlaps y. This  Ranges object will have additional metadata columns 
#' corresponding to the metadata columns in y. `find_overlaps_within()` is
#' the same but will only search for overlaps within y. For GRanges objects strand is
#' ignored, unless `find_overlaps_directed()` is used. If the Ranges objects have no
#' metadata, one could use `group_by_overlaps()` to be able to
#' identify the index of the input Range x that overlaps a Range in y. 
#' Alternatively,
#' `pair_overlaps()` could be used to place the x ranges next to the range
#' in y they overlap.
#'
#' @return A Ranges object with rows corresponding to the
#' ranges in x that overlap y.  In the case of `group_by_overlaps()`, returns
#' a GroupedRanges object, grouped by the number of overlaps
#' of ranges in x that overlap y (stored in a column called query).
#'
#' @examples
#' query <- data.frame(start = c(5,10, 15,20), width = 5, gc = runif(4)) %>%
#'              as_iranges()
#' subject <- data.frame(start = 2:6, width = 3:7, label = letters[1:5]) %>%
#'              as_iranges()
#'
#' find_overlaps(query, subject)
#' find_overlaps(query, subject, minoverlap = 5)
#' find_overlaps_within(query, subject) # same result as minoverlap
#' find_overlaps(query, subject, maxgap = 1)
#'
#' # -- GRanges objects, strand is ignored by default
#' query  <- data.frame(seqnames = "chr1",
#'                start = c(11,101),
#'                end = c(21, 200),
#'                name = c("a1", "a2"),
#'                strand = c("+", "-"),
#'                score = c(1,2)) %>%
#'            as_granges()
#' subject <- data.frame(seqnames = "chr1",
#'                       strand = c("+", "-", "+", "-"),
#'                       start = c(21,91,101,201),
#'                       end = c(30,101,110,210),
#'                       name = paste0("b", 1:4),
#'                       score = 1:4) %>%
#'                    as_granges()
#'
#' # ignores strandedness
#' find_overlaps(query, subject, suffix = c(".query", ".subject"))
#' find_overlaps(query, subject, suffix = c(".query", ".subject"), minoverlap = 2)
#' # adding directed prefix includes strand
#' find_overlaps_directed(query, subject, suffix = c(".query", ".subject"))
#'
#' @seealso [GenomicRanges::findOverlaps()], [IRanges::findOverlaps()]
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @export
find_overlaps <- function(x, y, maxgap, minoverlap, suffix = c(".x", ".y")) {
  UseMethod("find_overlaps")
}

#' @rdname ranges-overlaps
#' @export
find_overlaps.IntegerRanges <- function(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  .find_overlaps(x,y, suffix, findOverlaps, 
                 maxgap = maxgap, 
                 minoverlap = minoverlap, 
                 type = "any",
                 select = "all")
}

#' @rdname ranges-overlaps
#' @export
find_overlaps.GenomicRanges <- function(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  .find_overlaps(x,y, suffix, 
                 findOverlaps, 
                 maxgap = maxgap, 
                 minoverlap = minoverlap, 
                 type = "any",
                 select = "all",
                 ignore.strand = TRUE)
}

#' @rdname ranges-overlaps
#' @export
find_overlaps_within <- function(x, y, maxgap, minoverlap, suffix = c(".x", ".y")) {
  UseMethod("find_overlaps_within")
}

#' @rdname ranges-overlaps
#' @export
find_overlaps_within.IntegerRanges <- function(x,y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  .find_overlaps(x,y, suffix, findOverlaps, 
                 maxgap = maxgap, 
                 minoverlap = minoverlap, 
                 type = "within",
                 select = "all")
}

#' @rdname ranges-overlaps
#' @export
find_overlaps_within.GenomicRanges <- function(x,y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  .find_overlaps(x,y, suffix, findOverlaps, 
                 maxgap = maxgap, 
                 minoverlap = minoverlap, 
                 type = "within",
                 select = "all",
                 ignore.strand = TRUE)
}

#' @rdname ranges-overlaps
#' @export
find_overlaps_directed <- function(x, y, maxgap, minoverlap, suffix = c(".x", ".y")) {
  UseMethod("find_overlaps_directed")
}

#' @rdname ranges-overlaps
#' @export
find_overlaps_directed.GenomicRanges <- function(x,y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  .find_overlaps(x,y, suffix, findOverlaps, 
                 maxgap = maxgap, 
                 minoverlap = minoverlap, 
                 type = "any",
                 select = "all",
                 ignore.strand = FALSE)
}

#' @rdname ranges-overlaps
#' @export
find_overlaps_within_directed <- function(x, y, maxgap, minoverlap, suffix = c(".x", ".y")) {
  UseMethod("find_overlaps_within_directed")
}

#' @rdname ranges-overlaps
#' @export
find_overlaps_within_directed.GenomicRanges <- function(x, y, maxgap, minoverlap, suffix = c(".x", ".y")) {
  .find_overlaps(x,y, suffix, findOverlaps, 
                 maxgap = maxgap, 
                 minoverlap = minoverlap, 
                 type = "within",
                 select = "all",
                 ignore.strand = FALSE)
}
