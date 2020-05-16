#' Find nearest neighbours between two Ranges objects
#'
#' @param x,y Ranges objects, add the nearest neighbours of ranges in x to
#' those in y.
#' @param suffix A character vector of length two used to identify metadata columns
#' @param distance logical vector whether to add "distance" column to output. If set to
#'   a character vector of length 1, will use that as distance column name. If a
#'   column matching the distance column already exits, will prepend the suffix
#'   for y with a warning.
#'
#' @details By default `join_nearest` will find arbitrary nearest
#' neighbours in either direction and ignore any strand information.
#' The `join_nearest_left` and `join_nearest_right`  methods
#' will find arbitrary nearest neighbour ranges on x that are left/right of
#' those on y and ignore any strand information.
#'
#' The `join_nearest_upstream` method will find arbitrary nearest
#' neighbour ranges on x that are upstream of those on y. This takes into
#' account strandedness of the ranges.
#' On the positive strand nearest upstream will be on the
#' left and on the negative strand nearest upstream will be on the right.
#'
#' The `join_nearest_downstream` method will find arbitrary nearest
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
join_nearest <- function(x, y, suffix = c(".x", ".y"), distance = FALSE) {
  UseMethod("join_nearest")
}

#' @export
join_nearest.IntegerRanges <- function(x,y, suffix = c(".x", ".y"), distance = FALSE) {
  hits <- make_hits(x, y, distanceToNearest, select = "arbitrary")
  join <- expand_by_hits(x, y, suffix, hits)
  add_nearest_hits_distance(join, hits, suffix[2], distance)
}

#' @export
join_nearest.GenomicRanges <- function(x,y, suffix = c(".x", ".y"), distance = FALSE) {
  hits <- make_hits(x, y, distanceToNearest, select = "arbitrary", ignore.strand = TRUE)
  join <- expand_by_hits(x, y, suffix, hits)
  add_nearest_hits_distance(join, hits, suffix[2], distance)
}

#' @rdname ranges-nearest
#' @export
join_nearest_left <- function(x, y, suffix = c(".x", ".y"), distance = FALSE) {
  UseMethod("join_nearest_left")
}

#' @export
join_nearest_left.IntegerRanges <- function(x,y, suffix = c(".x", ".y"), distance = FALSE) {
  hits <- make_hits(x, y, distanceToNearest, select = "all")
  mcols(hits)$is_left <- start(y[subjectHits(hits)]) <= start(x[queryHits(hits)]) &
    end(y[subjectHits(hits)]) <= start(x[queryHits(hits)])
  hits <- hits[mcols(hits)$is_left]
  join <- expand_by_hits(x,y, suffix, hits)
  add_nearest_hits_distance(join, hits, suffix[2], distance)
}

#' @export
join_nearest_left.GenomicRanges <- function(x,y, suffix = c(".x", ".y"), distance = FALSE) {
  hits <- make_hits(x, y, distanceToNearest, select = "all", ignore.strand = TRUE)
  mcols(hits)$is_left <- start(y[subjectHits(hits)]) <= start(x[queryHits(hits)]) &
    end(y[subjectHits(hits)]) <= start(x[queryHits(hits)])
  hits <- hits[mcols(hits)$is_left]
  join <- expand_by_hits(x,y, suffix, hits)
  add_nearest_hits_distance(join, hits, suffix[2], distance)
}

#' @importFrom IRanges nearest
#' @rdname ranges-nearest
#' @export
join_nearest_right <- function(x, y,  suffix = c(".x", ".y"), distance = FALSE) { UseMethod("join_nearest_right")}

#' @export
join_nearest_right.IntegerRanges <- function(x, y, suffix = c(".x", ".y"), distance = FALSE) {
  hits <- make_hits(x, y, distanceToNearest, select = "all")
  mcols(hits)$is_right <- end(x[queryHits(hits)]) <= start(y[subjectHits(hits)])
  hits <- hits[mcols(hits)$is_right]
  join <- expand_by_hits(x,y, suffix, hits)
  add_nearest_hits_distance(join, hits, suffix[2], distance)
}

#' @export
join_nearest_right.GenomicRanges <- function(x, y,  suffix = c(".x", ".y"), distance = FALSE) {
  hits <- make_hits(x, y, distanceToNearest, select = "all", ignore.strand = TRUE)
  mcols(hits)$is_right <- end(x[queryHits(hits)]) <= start(y[subjectHits(hits)])
  hits <- hits[mcols(hits)$is_right]
  join <- expand_by_hits(x,y, suffix, hits)
  add_nearest_hits_distance(join, hits, suffix[2], distance)
}


#' @rdname ranges-nearest
#' @export
join_nearest_upstream <- function(x, y,  suffix = c(".x", ".y"), distance = FALSE) { UseMethod("join_nearest_upstream")}

#' @export
join_nearest_upstream.GenomicRanges <- function(x, y,  suffix = c(".x", ".y"), distance = FALSE) {
  hits <- distanceToNearest(x,y, select = "all", ignore.strand = FALSE)
  mcols(hits)$is_right <- end(x[queryHits(hits)]) <= start(y[subjectHits(hits)])
  mcols(hits)$is_left <- start(x[queryHits(hits)]) >= start(y[subjectHits(hits)])
  mcols(hits)$direction <- strand(x[queryHits(hits)])
  mcols(hits)$is_upstream <- ifelse(mcols(hits)$direction == "+",
                                    mcols(hits)$is_left,
                                    mcols(hits)$is_right)
  hits <- hits[mcols(hits)$is_upstream]
  join <- expand_by_hits(x,y, suffix, hits)
  add_nearest_hits_distance(join, hits, suffix[2], distance)
}

#' @rdname ranges-nearest
#' @importFrom IRanges nearest
#' @export
join_nearest_downstream <- function(x, y,  suffix = c(".x", ".y"), distance = FALSE) { UseMethod("join_nearest_downstream")}

#' @export
join_nearest_downstream.GenomicRanges <- function(x, y, suffix = c(".x", ".y"), distance = FALSE) {
  hits <- distanceToNearest(x,y, select = "all", ignore.strand = FALSE)
  mcols(hits)$is_right <- end(x[queryHits(hits)]) <= start(y[subjectHits(hits)])
  mcols(hits)$is_left <- start(x[queryHits(hits)]) >= start(y[subjectHits(hits)])
  mcols(hits)$direction <- strand(x[queryHits(hits)])
  # on positive strand nearest downtream will be on the right
  # on negative strand nearest downstream will be on the left
  mcols(hits)$is_downstream <- ifelse(mcols(hits)$direction == "-",
                                    mcols(hits)$is_left,
                                    mcols(hits)$is_right)

  hits <- hits[mcols(hits)$is_downstream]
  join <- expand_by_hits(x,y, suffix, hits)
  add_nearest_hits_distance(join, hits, suffix[2], distance)
}

# -----------------------
# Nearest Distance Utils

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

#' Add suffix to column name
#' 
#' I think dplyr::add_count will change how they add names in 1.0.0
#' used to it would add n until it got a unique name
#' I can't find the blogpost describing the new behavior,
#' but it may be good to mirror whatever that is here.
#' For now, just recursively add suffix until a unique column is found
#'
#' @param name column name to add
#' @param suffix suffix to use
#' @param names names of data.frame in which to add new column
#'
#' @return a unique name generated by prepending suffix recursively until
#'     a unique name is found
#' @noRd
add_suffix <- function(name, suffix, names){
  
  stopifnot(length(suffix) == 1)
  
  if (name %in% names) {
    name <- paste0(name, suffix)
  } else {
    return(name)
  }
  add_suffix(name, suffix, names)
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

#' @param ranges ranges object from join_nearest expand_hits internal function
#'
#' @param hits  hits from distanceToNearest
#' @param colname name of column to hold distance values
#' @param suffix suffix to add to colname if exists
#'
#' @noRd
add_distance_col <- function(ranges, hits, colname = "distance", suffix = ".y"){
  UseMethod("add_distance_col")
}

#' @noRd
add_distance_col.GenomicRanges <- function(ranges, hits, colname = "distance", suffix = ".y"){
  
  if (!("distance" %in% names(mcols(hits)))) {
    stop("hits object must contain a distance column")
  }
  
  colname <- add_suffix(colname, suffix, names(mcols(ranges)))
  
  mcols(ranges)[colname] <- mcols(hits)$distance
  
  return(ranges)
}

#' @noRd
add_distance_col.IRanges <- function(ranges, hits, colname = "distance", suffix = ".y"){
  
  if (!("distance" %in% names(mcols(hits)))) {
    stop("hits object must contain a distance column")
  }
  
  colname <- add_suffix(colname, suffix = suffix, names(mcols(ranges)))
  
  if (is.null(mcols(ranges))){
    metadata <- DataFrame(colname = mcols(hits)$distance)
    mcols(ranges) <- metadata
  } else {
    mcols(ranges)[colname] <- mcols(hits)$distance
  }
  return(ranges)
}