#' Find nearest neighbours between two Ranges objects
#'
#' @param x,y Ranges objects, add the nearest neighbours of ranges in x to
#' those in y.
#' @param suffix A character vector of length two used to identify metadata columns
#' coming from x and y.
#'
#' @details By default \code{join_nearest} will find abritrary nearest
#' neighbours in either direction and ignore any strand information.
#' The \code{join_nearest_left} and \code{join_nearest_right}  methods
#' will find abritrary nearest neighbour ranges on x that are left/right of
#' those on y and ignore any strand information.
#'
#' The \code{join_nearest_upstream} method will find abritrary nearest
#' neighbour ranges on x that are upstream of those on y. This takes into
#' account strandedness of the ranges.
#' On the positive strand nearest upstream will be on the
#' left and on the negative strand nearest upstream will be on the right.
#'
#' The \code{join_nearest_downstream} method will find abritrary nearest
#' neighbour ranges on x that are upstream of those on y. This takes into
#' account strandedness of the ranges.On the positive strand nearest downstream
#' will be on the right and on the negative strand nearest upstream will be on
#' the left.
#'
#' @return A Ranges object corresponding to the nearest ranges, all metadata
#' is copied over from the right-hand side ranges `y`.
#'
#' @examples
#' query <- data.frame(start = c(5,10, 15,20),
#'                    width = 5,
#'                    gc = runif(4)) %>%
#'              as_iranges()
#' subject <- data.frame(start = c(2:6, 24),
#'                       width = 3:8,
#'                       label = letters[1:6]) %>%
#'              as_iranges()
#'
#' join_nearest(query, subject)
#' join_nearest_left(query, subject)
#' join_nearest_right(query, subject)
#'
#' subject  <- data.frame(seqnames = "chr1",
#'                start = c(11,101),
#'                end = c(21, 200),
#'                name = c("a1", "a2"),
#'                strand = c("+", "-"),
#'                score = c(1,2)) %>%
#'            as_granges()
#' query <- data.frame(seqnames = "chr1",
#'                       strand = c("+", "-", "+", "-"),
#'                       start = c(21,91,101,201),
#'                       end = c(30,101,110,210),
#'                       name = paste0("b", 1:4),
#'                       score = 1:4) %>%
#'                    as_granges()
#' join_nearest_upstream(query, subject)
#' join_nearest_downstream(query, subject)
#' @rdname ranges-nearest
#' @importFrom IRanges nearest
#' @export
join_nearest <- function(x, y, suffix = c(".x", ".y")) {
  UseMethod("join_nearest")
}

#' @export
join_nearest.Ranges <- function(x,y, suffix = c(".x", ".y")) {
  hits <- nearest(x,y, select = "arbitrary")
  no_hits_id <- !is.na(hits)
  left <- x[no_hits_id, ]
  right <- y[hits[no_hits_id], ]
  mcols(left) <- mcols_overlaps_update(left, right, suffix)
  return(left)
}

#' @export
join_nearest.GenomicRanges <- function(x,y, suffix = c(".x", ".y")) {
  hits <- nearest(x,y, select = "arbitrary", ignore.strand = TRUE)
  no_hits_id <- !is.na(hits)
  left <- x[no_hits_id, ]
  right <- y[hits[no_hits_id], ]
  mcols(left) <- mcols_overlaps_update(left, right, suffix)
  return(left)
}

#' @rdname ranges-nearest
#' @export
join_nearest_left <- function(x, y, suffix = c(".x", ".y")) {
  UseMethod("join_nearest_left")
}

#' @export
join_nearest_left.Ranges <- function(x,y, suffix = c(".x", ".y")) {
  hits <- nearest(x,y, select = "all")
  mcols(hits)$is_left <- start(y[subjectHits(hits)]) <= start(x[queryHits(hits)]) &
    end(y[subjectHits(hits)]) <= start(x[queryHits(hits)])
  hits <- hits[mcols(hits)$is_left]
  left <- x[queryHits(hits), ]
  right <- y[subjectHits(hits), ]
  mcols(left) <- mcols_overlaps_update(left, right, suffix)
  left
}

#' @export
join_nearest_left.GenomicRanges <- function(x,y, suffix = c(".x", ".y")) {
  hits <- nearest(x,y, select = "all", ignore.strand = TRUE)
  mcols(hits)$is_left <- start(y[subjectHits(hits)]) <= start(x[queryHits(hits)]) &
    end(y[subjectHits(hits)]) <= start(x[queryHits(hits)])
  hits <- hits[mcols(hits)$is_left]
  left <- x[queryHits(hits), ]
  right <- y[subjectHits(hits), ]
  mcols(left) <- mcols_overlaps_update(left, right, suffix)
  left
}

#' @importFrom IRanges nearest
#' @rdname ranges-nearest
#' @export
join_nearest_right <- function(x, y,  suffix = c(".x", ".y")) { UseMethod("join_nearest_right")}

#' @export
join_nearest_right.Ranges <- function(x, y, suffix = c(".x", ".y")) {
  hits <- nearest(x,y, select = "all")
  mcols(hits)$is_right <- end(x[queryHits(hits)]) <= start(y[subjectHits(hits)])
  hits <- hits[mcols(hits)$is_right]
  left <- x[queryHits(hits), ]
  right <- y[subjectHits(hits), ]
  mcols(left) <- mcols_overlaps_update(left, right, suffix)
  left
}

#' @export
join_nearest_right.GenomicRanges <- function(x, y,  suffix = c(".x", ".y")) {
  hits <- nearest(x,y, select = "all", ignore.strand = TRUE)
  mcols(hits)$is_right <- end(x[queryHits(hits)]) <= start(y[subjectHits(hits)])
  hits <- hits[mcols(hits)$is_right]
  left <- x[queryHits(hits), ]
  right <- y[subjectHits(hits), ]
  mcols(left) <- mcols_overlaps_update(left, right, suffix)
  left
}


#' @rdname ranges-nearest
#' @export
join_nearest_upstream <- function(x, y,  suffix = c(".x", ".y")) { UseMethod("join_nearest_upstream")}

#' @export
join_nearest_upstream.GenomicRanges <- function(x, y,  suffix = c(".x", ".y")) {
  hits <- nearest(x,y, select = "all", ignore.strand = FALSE)
  mcols(hits)$is_right <- end(x[queryHits(hits)]) <= start(y[subjectHits(hits)])
  mcols(hits)$is_left <- start(x[queryHits(hits)]) >= start(y[subjectHits(hits)])
  mcols(hits)$direction <- strand(x[queryHits(hits)])
  mcols(hits)$is_upstream <- ifelse(mcols(hits)$direction == "+",
                                    mcols(hits)$is_left,
                                    mcols(hits)$is_right)
  hits <- hits[mcols(hits)$is_upstream]
  left <- x[queryHits(hits), ]
  right <- y[subjectHits(hits), ]
  mcols(left) <- mcols_overlaps_update(left, right, suffix)
  left
}

#' @rdname ranges-nearest
#' @importFrom IRanges nearest
#' @export
join_nearest_downstream <- function(x, y,  suffix = c(".x", ".y")) { UseMethod("join_nearest_downstream")}

#' @export
join_nearest_downstream.GenomicRanges <- function(x, y, suffix = c(".x", ".y")) {
  hits <- nearest(x,y, select = "all", ignore.strand = FALSE)
  mcols(hits)$is_right <- end(x[queryHits(hits)]) <= start(y[subjectHits(hits)])
  mcols(hits)$is_left <- start(x[queryHits(hits)]) >= start(y[subjectHits(hits)])
  mcols(hits)$direction <- strand(x[queryHits(hits)])
  # on positive strand nearest downtream will be on the right
  # on negative strand nearest downstream will be on the left
  mcols(hits)$is_downstream <- ifelse(mcols(hits)$direction == "-",
                                    mcols(hits)$is_left,
                                    mcols(hits)$is_right)

  hits <- hits[mcols(hits)$is_downstream]
  left <- x[queryHits(hits), ]
  right <- y[subjectHits(hits), ]
  mcols(left) <- mcols_overlaps_update(left, right, suffix)
  left
}
