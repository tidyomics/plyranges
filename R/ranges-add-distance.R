#' @param ranges ranges object from join_nearest expand_hits internal function
#' @param hits  hits from distanceToNearest
#' @param colname name of column to hold distance values
#'
#' @noRd
add_distance_col <- function(ranges, hits, colname) {
  
  if (colname %in% names(mcols(ranges))){
    stop(paste0(colname, " already exists in destination metadata"))
  }
  
  if (is.null(mcols(ranges))){
    # handle IRanges NULL adding X column of NA's
    meta <- DataFrame("distance" = NA_integer_)
    names(meta) <- colname
    mcols(ranges) <- meta
  } else {
    mcols(ranges)[[colname]] <-  NA_integer_
  }
  
  mcols(ranges[queryHits(hits)])[[colname]] <- mcols(hits)$distance
  
  return(ranges)
}

#' Macro for building add_nearest_distance_* family of functions
#' Change this to edit defaults for all functions
#' @param fun nearest hits helper function
#' @noRd
make_add_nearest_distance <- function(fun){
  f <- function(x, y = x, name = "distance"){
    hits <- fun(x, y)
    
    add_distance_col(x, hits, colname = name)
  }
}

#' Add distance to nearest neighbours between two Ranges objects
#' 
#' Appends distance to nearest subject range to query ranges similar to setting
#' `distance` in `join_nearest_`. Distance is set to `NA` for features with no
#' nearest feature by the selected nearest metric.
#'
#' @param x The query ranges
#' @param y the subject ranges within which the nearest ranges are found.
#'   If missing, query ranges are used as the subject.
#' @param name column name to create containing distance values
#'   
#' @details By default `add_nearest_distance` will find arbitrary nearest
#' neighbours in either direction and ignore any strand information.
#' The `add_nearest_distance_left` and `add_nearest_distance_right`  methods
#' will find arbitrary nearest neighbour ranges on x that are left/right of
#' those on y and ignore any strand information.
#'
#' The `add_nearest_distance_upstream` method will find arbitrary nearest
#' neighbour ranges on x that are upstream of those on y. This takes into
#' account strandedness of the ranges.
#' On the positive strand nearest upstream will be on the
#' left and on the negative strand nearest upstream will be on the right.
#'
#' The `add_nearest_distance_downstream` method will find arbitrary nearest
#' neighbour ranges on x that are upstream of those on y. This takes into
#' account strandedness of the ranges. On the positive strand nearest downstream
#' will be on the right and on the negative strand nearest upstream will be on
#' the left.
#'
#' @return ranges in `x` with additional column containing the distance to the
#'   nearest range in `y`.
#'  
#' @rdname add-nearest-distance
#' @seealso \code{\link{join_nearest}}
#' @export
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
#' add_nearest_distance(query, subject)
#' add_nearest_distance_left(query, subject)
#' add_nearest_distance_left(query)
#' @export
add_nearest_distance <- make_add_nearest_distance(hits_nearest)

#' @rdname add-nearest-distance
#' @export
add_nearest_distance_left <- make_add_nearest_distance(hits_nearest_left)

#' @rdname add-nearest-distance
#' @export
add_nearest_distance_right <- make_add_nearest_distance(hits_nearest_right)

#' @rdname add-nearest-distance
#' @export
add_nearest_distance_upstream <- make_add_nearest_distance(hits_nearest_upstream)

#' @rdname add-nearest-distance
#' @export
add_nearest_distance_downstream <- make_add_nearest_distance(hits_nearest_downstream)

