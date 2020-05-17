#' Macro for building add_nearest_distance_* family of functions
#' Change this to edit defaults for all functions
make_add_nearest_distance <- function(fun){
  f <- function(query, subject = query, .id = "distance", suffix = ".y"){
    hits <- fun(query, subject)
    
    add_distance_col(query, hits, colname = .id, suffix = suffix)
  }
}

#' Add distance to nearest neighbours between two Ranges objects
#' 
#' Appends distance to nearest subject range to query ranges similar to setting
#' `distance` in `join_nearest_`. Distance is set to `NA` for features with no
#' nearest feature by the selected nearest metric.
#'
#' @param query The query ranges
#' @param subject the subject ranges within which the nearest ranges are found.
#'   If missing, query ranges are used as the subject.
#' @param .id column name to create containing distance values
#' @param suffix if .id already exists as column in query, prepend this value to .id column name
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
#' @return query ranges with additional column containing the distance to the
#'   nearest range in subject.
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

#' @param expanded_hits output of expand_hits()
#' @param hits output from make_hits() when f_in() = nearestDistance
#' @param suffix suffix to prepend to distance column if name exists already
#' @param distance `logical(1)` whether to add "distance" column to output. If set to
#'   `character(1)`, will use that as distance column name. (Default: FALSE)
#'
#' @noRd
add_nearest_hits_distance <- function(expanded_hits, hits, suffix = ".y", distance = FALSE){
  
  if (length(distance) > 1){
    stop("distance must be of length 1")
  }
  
  if (distance == TRUE){
    
    expanded_hits <- add_distance_col(expanded_hits, hits, colname = "distance", suffix = suffix)
    
  } else if (is.character(distance)) {
    
    expanded_hits <- add_distance_col(expanded_hits, hits, colname = distance, suffix = suffix)
  }
  
  return(expanded_hits)
}