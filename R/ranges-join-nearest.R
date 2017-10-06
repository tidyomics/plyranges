bind_nearest <- function(x, y, type = "nearest") {

  if (ncol(y) > 1) {
    names(y) <- c(type, paste0(names(y)[-1], ".y"))
  } else {
    names(y) <- type
  }

  if (is.null(mcols(x)) & is(x, "IRanges")) {
    mcols(x) <- y
  } else {
    mcols(x) <- cbind(mcols(x), y)
  }
  return(x)
}

nearest_rng <- function(x, y, hits, type = "nearest") {
  no_hits_id <- !is.na(hits)
  expand_y <- as(y[hits[no_hits_id]], "DataFrame")
  reduce_x <- x[no_hits_id]
  bind_nearest(reduce_x, expand_y, type)
}

nearest_rng_all <- function(x, y, hits, type = "nearest") {
  reduce_x <- x[queryHits(hits)]
  expand_y <- as(y[subjectHits(hits)], "DataFrame")
  bind_nearest(reduce_x, expand_y, type)
}
#' Find nearest neighbours between two Ranges objects
#'
#' @param x,y Ranges objects, add the nearest neighbours of ranges in x to
#' those in y.
#'
#' @details By default \code{join_nearest} will find abritrary nearest
#' neighbours in either direction and ignore any strand information.
#'
#' @return A Ranges object with a metadata column called nearest that
#' contains the corresponding nearest Ranges.
#' @rdname nearest-ranges
#' @importFrom IRanges nearest
#' @export
join_nearest <- function(x, y) {
  UseMethod("join_nearest")
}

#' @export
join_nearest.Ranges <- function(x,y) {
  hits <- nearest(x,y, select = "arbitrary")
  nearest_rng(x,y, hits)
}

#' @export
join_nearest.GenomicRanges <- function(x,y) {
  hits <- nearest(x,y, select = "arbitrary", ignore.strand = FALSE)
  nearest_rng(x,y, hits)
}

#' Find nearest neighbours on the left between two Ranges objects
#'
#' @param x,y Ranges objects, add the nearest neighbours of ranges in x that
#' are left of those in y.
#'
#' @details By default \code{join_nearest_left} will find abritrary nearest
#' neighbour ranges on x that are left of those on y and ignore any
#' strand information.
#'
#' @return A Ranges object with a metadata column called nearest that
#' contains the corresponding nearest Ranges.
#' @rdname nearest-left-ranges
#' @importFrom IRanges nearest
#' @export
#' @export
join_nearest_left <- function(x, y) { UseMethod("join_nearest_left")}

#' @export
join_nearest_left.Ranges <- function(x,y) {
  hits <- nearest(x,y, select = "all")
  mcols(hits)$is_left <- start(x[queryHits(hits)]) >= start(y[subjectHits(hits)]) &
    end(x[queryHits(hits)]) > start(y[subjectHits(hits)])
  hits <- hits[mcols(hits)$is_left]
  nearest_rng_all(x,y, hits)
}

#' @export
join_nearest_left.GenomicRanges <- function(x,y) {
  hits <- nearest(x,y, select = "all", ignore.strand = TRUE)
  mcols(hits)$is_left <- start(x[queryHits(hits)]) >= start(y[subjectHits(hits)])
  hits <- hits[mcols(hits)$is_left]
  nearest_rng_all(x,y, hits)

}

#' Find nearest neighbours on the right between two Ranges objects
#'
#' @param x,y Ranges objects, add the nearest neighbours of ranges in x that
#' are right of those in y.
#'
#' @details By default \code{join_nearest_right} will find abritrary nearest
#' neighbour ranges on x that are right of those on y and ignore any
#' strand information.
#'
#' @return A Ranges object with a metadata column called nearest that
#' contains the corresponding nearest Ranges.
#' @rdname nearest-right-ranges
#' @importFrom IRanges nearest
#' @export
#' @export
join_nearest_right <- function(x, y) { UseMethod("join_nearest_right")}

#' @export
join_nearest_right.Ranges <- function(x, y) {
  hits <- nearest(x,y, select = "all")
  mcols(hits)$is_right <- end(x[queryHits(hits)]) <= start(y[subjectHits(hits)])
  hits <- hits[mcols(hits)$is_right]
  nearest_rng_all(x, y, hits)

}

#' @export
join_nearest_right.GenomicRanges <- function(x, y) {
  hits <- nearest(x,y, select = "all", ignore.strand = TRUE)
  mcols(hits)$is_right <- end(x[queryHits(hits)]) <= start(y[subjectHits(hits)])
  hits <- hits[mcols(hits)$is_right]
  nearest_rng_all(x,y, hits)
}

#' Find nearest neighbours that are upstream between two Ranges objects
#'
#' @param x,y GRanges objects, add the nearest neighbours of ranges in x that
#' are upstream of those in y.
#'
#' @details \code{join_nearest_upstream} will find abritrary nearest
#' neighbour ranges on x that are upstream of those on y. This takes into
#' account strandedness of the ranges.
#' On the positive strand nearest upstream will be on the
#' left and on the negative strand nearest upstream will be on the right.
#'
#' @return A Ranges object with a metadata column called nearest that
#' contains the corresponding nearest Ranges.
#' @rdname nearest-upstream-ranges
#' @importFrom IRanges nearest
#' @export
join_nearest_upstream <- function(x, y) { UseMethod("join_nearest_upstream")}

#' @export
join_nearest_upstream.GenomicRanges <- function(x, y) {
  hits <- nearest(x,y, select = "all", ignore.strand = FALSE)
  mcols(hits)$is_right <- end(x[queryHits(hits)]) <= start(y[subjectHits(hits)])
  mcols(hits)$is_left <- start(x[queryHits(hits)]) >= start(y[subjectHits(hits)])
  mcols(hits)$direction <- strand(x[queryHits(hits)])
  mcols(hits)$is_upstream <- ifelse(mcols(hits)$direction == "+",
                                    mcols(hits)$is_left,
                                    mcols(hits)$is_right)
  hits <- hits[mcols(hits)$is_upstream]
  nearest_rng_all(x,y, hits)

}


#' Find nearest neighbours that are downstream between two Ranges objects
#'
#' @param x,y GRanges objects, add the nearest neighbours of ranges in x that
#' are downstream of those in y.
#'
#' @details \code{join_nearest_downstream} will find abritrary nearest
#' neighbour ranges on x that are upstream of those on y. This takes into
#' account strandedness of the ranges.On the positive strand nearest downstream
#' will be on the right and on the negative strand nearest upstream will be on
#' the left.
#'
#' @return A Ranges object with a metadata column called nearest that
#' contains the corresponding nearest Ranges.
#' @rdname nearest-downstream-ranges
#' @importFrom IRanges nearest
#' @export
join_nearest_downstream <- function(x, y) { UseMethod("join_nearest_downstream")}

#' @export
join_nearest_downstream.GenomicRanges <- function(x, y) {
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
  nearest_rng_all(x,y,hits)
}
