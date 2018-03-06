#' Pair together two ranges objects
#'
#' @param x,y Ranges objects to pair together.
#' @param maxgap,minoverlap The maximimum gap between intervals as an integer
#' greater than or equal to negative one. The minimum amount of overlap between intervals
#' as an integer greater than zero, accounting for the maximum gap.
#' @param suffix A character vector of length two used to identify metadata columns
#' coming from x and y.
#'
#' @details These functions return a DataFrame object, and is one way of
#' representing paired alignments with plyranges.
#'
#' @return a DataFrame with two ranges columns and the corresponding metadata
#' columns.
#'
#' @seealso [join_nearest()][join_overlap_inner()][join_precede()][join_follow()]
#' @examples
#' #' query <- data.frame(start = c(5,10, 15,20), width = 5, gc = runif(4)) %>%
#'              as_iranges()
#' subject <- data.frame(start = 2:6, width = 3:7, label = letters[1:5]) %>%
#'              as_iranges()
#'
#' pair_overlaps(query, subject)
#' pair_overlaps(query, subject, minoverlap = 5)
#' pair_nearest(query, subject)
#'
#'
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
#' pair_overlaps(query, subject, suffix = c(".query", ".subject"))
#' pair_follow(query, subject, suffix = c(".query", ".subject"))
#' pair_precede(query, subject, suffix = c(".query", ".subject"))
#' pair_precede(query, subject, suffix = c(".query", ".subject"))
#' @rdname ranges-pairs
#' @export
pair_overlaps <- function(x, y, maxgap, minoverlap, suffix) {
  UseMethod("pair_overlaps")
}

#' @export
pair_overlaps.Ranges <- function(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  hits <- findOverlaps(x,y, maxgap, minoverlap, type = "any", select = "all")
  left <- x[queryHits(hits), ]
  right <- y[subjectHits(hits), ]
  mcols_overlaps_update(left, right, suffix, return_data_frame = TRUE)
}

#' @export
pair_overlaps.GenomicRanges <- function(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  hits <- findOverlaps(x,y, maxgap, minoverlap,
                       type = "any", select = "all", ignore.strand = TRUE)
  left <- x[queryHits(hits), ]
  right <- y[subjectHits(hits), ]
  mcols_overlaps_update(left, right, suffix, return_data_frame = TRUE)
}

#' @rdname ranges-pairs
#' @export
pair_nearest <- function(x, y, suffix) {
  UseMethod("pair_nearest")
}

#' @export
pair_nearest.Ranges <- function(x,y, suffix = c(".x", ".y")) {
  hits <- nearest(x,y, select = "arbitrary")
  no_hits_id <- !is.na(hits)
  left <- x[no_hits_id, ]
  right <- y[hits[no_hits_id], ]
  mcols_overlaps_update(left, right, suffix, return_data_frame = TRUE)
}

#' @export
pair_nearest.GenomicRanges <- function(x,y, suffix = c(".x", ".y")) {
  hits <- nearest(x,y, select = "arbitrary", ignore.strand = TRUE)
  no_hits_id <- !is.na(hits)
  left <- x[no_hits_id, ]
  right <- y[hits[no_hits_id], ]
  mcols_overlaps_update(left, right, suffix, return_data_frame = TRUE)
}

#' @rdname ranges-pairs
#' @export
pair_precede <- function(x, y, suffix) {
  UseMethod("pair_precede")
}

#' @export
pair_precede.Ranges <- function(x,y, suffix = c(".x", ".y")) {
  hits <- precede(x,y)
  no_hits_id <- !is.na(hits)
  left <- x[no_hits_id, ]
  right <- y[hits[no_hits_id], ]
  mcols_overlaps_update(left, right, suffix, return_data_frame = TRUE)
}

#' @export
pair_precede.GenomicRanges <- function(x,y, suffix = c(".x", ".y")) {
  hits <- precede(x,y, ignore.strand = TRUE)
  no_hits_id <- !is.na(hits)
  left <- x[no_hits_id, ]
  right <- y[hits[no_hits_id], ]
  mcols_overlaps_update(left, right, suffix, return_data_frame = TRUE)
}

#' @rdname ranges-pairs
#' @export
pair_follow <- function(x, y, suffix) {
  UseMethod("pair_follow")
}
#' @export
pair_follow.Ranges <- function(x,y, suffix = c(".x", ".y")) {
  hits <- follow(x,y)
  no_hits_id <- !is.na(hits)
  left <- x[no_hits_id, ]
  right <- y[hits[no_hits_id], ]
  mcols_overlaps_update(left, right, suffix, return_data_frame = TRUE)
}

#' @export
pair_follow.GenomicRanges <- function(x,y, suffix = c(".x", ".y")) {
  hits <- follow(x,y, ignore.strand = TRUE)
  no_hits_id <- !is.na(hits)
  left <- x[no_hits_id, ]
  right <- y[hits[no_hits_id], ]
  mcols_overlaps_update(left, right, suffix, return_data_frame = TRUE)
}
