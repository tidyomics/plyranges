#' Add distance to nearest feature
#' 
#' A light wrapper around `distanceToNearest()` for appending distance values to query ranges.
#'
#' @param query The query ranges
#' @param subject the subject ranges within which the nearest ranges are found.
#'   If missing, query ranges are used as the subject.
#' @param ... Additional arguments passed to distanceToNearest()
#' @param .id column name to create containing distance values
#' @param suffix if .id already exists as column in query, prepend this value to .id column name
#'
#' @return query ranges with additional column containing the distance to the
#'   nearest range in subject.
#'  
#' @seealso ranges-nearest
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
#' add_nearest_distance(query, subject, ignore.strand = TRUE)
add_nearest_distance <- function(query, subject = query, ..., .id = "distance", suffix = ".y"){
  
  hits <- distanceToNearest(query, subject, ...)
  
  add_distance_col(query, hits, colname = .id, suffix = suffix)
  
}

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