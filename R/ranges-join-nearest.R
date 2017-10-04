
join_nearest <- function(x, y) {
  UseMethod("join_nearest")
}

join_nearest.Ranges <- function(x,y) {
  hits <- nearest(x,y, select = "arbitrary")
  no_hits_id <- !is.na(hits)
  expand_y <- as(y[hits[no_hits_id]], "DataFrame")

  if (ncol(expand_y) > 1) {
    names(expand_y) <- c("nearest", paste0(names(expand_y)[-1], ".y"))
  } else {
    names(expand_y) <- c("nearest")
  }

  reduce_x <- x[no_hits_id]
  if (is.null(mcols(reduce_x))) {
    mcols(reduce_x) <- expand_y
  } else {
    mcols(reduce_x) <- cbind(mcols(reduce_x), expand_y)
  }

  return(reduce_x)
}

join_nearest.GenomicRanges <- function(x,y) {
  hits <- nearest(x,y, select = "arbitrary", ignore.strand = FALSE)
  no_hits_id <- !is.na(hits)
  expand_y <- as(y[hits[no_hits_id]], "DataFrame")
  if (ncol(expand_y) > 1) {
    names(expand_y) <- c("nearest", paste0(names(expand_y)[-1], ".y"))
  } else {
    names(expand_y) <- c("nearest")
  }
  reduce_x <- x[no_hits_id]
  mcols(reduce_x) <- cbind(mcols(reduce_x), expand_y)
  return(reduce_x)
}

join_nearest_left <- function(x, y) { UseMethod("join_nearest_left")}

join_nearest_left.Ranges <- function(x,y) {
  hits <- nearest(x,y, select = "all")
  mcols(hits)$is_left <- start(x[queryHits(hits)]) <= start(y[subjectHits(hits)]) &
    end(x[queryHits(hits)]) > start(y[subjectHits(hits)])
  #mcols(hits)$dist <- distance(x[queryHits(hits)], y[subjectHits(hits)])
  hits <- hits[mcols(hits)$is_left]
  # should we keep all hits here or select arbitary?
  reduce_x <- x[queryHits(hits)]
  expand_y <- as(y[subjectHits(hits)], "DataFrame")
  if (ncol(expand_y) > 1) {
    names(expand_y) <- c("nearest", paste0(names(expand_y)[-1], ".y"))
  } else {
    names(expand_y) <- c("nearest")
  }

  if (is.null(mcols(reduce_x))) {
    mcols(reduce_x) <- expand_y
  } else {
    mcols(reduce_x) <- cbind(mcols(reduce_x), expand_y)
  }
  reduce_x
}

join_nearest_right <- function(x, y) { UseMethod("join_nearest_right")}

join_nearest_right.Ranges <- function(x, y) {
  hits <- nearest(x,y, select = "all")
  mcols(hits)$is_right <- end(x[queryHits(hits)]) <= start(y[subjectHits(hits)])
  hits <- hits[mcols(hits)$is_right]
  # should we keep all hits here or select arbitary?
  reduce_x <- x[queryHits(hits)]
  expand_y <- as(y[subjectHits(hits)], "DataFrame")
  if (ncol(expand_y) > 1) {
    names(expand_y) <- c("nearest", paste0(names(expand_y)[-1], ".y"))
  } else {
    names(expand_y) <- c("nearest")
  }

  if (is.null(mcols(reduce_x))) {
    mcols(reduce_x) <- expand_y
  } else {
    mcols(reduce_x) <- cbind(mcols(reduce_x), expand_y)
  }
  reduce_x
}

join_nearest_upstream <- function(x, y) {

}

join_nearest_downstream <- function(x, y) {

}
