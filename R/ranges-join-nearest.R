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
join_nearest <- function(x,y, suffix = c(".x", ".y"), distance = FALSE) {
  hits <- hits_nearest(x, y)
  distance_col <- set_distance_col(distance)
  expand_by_hits(x, y, suffix, hits, hits_mcols_to_keep = distance_col)
}

#' @rdname ranges-nearest
#' @export
join_nearest_left <- function(x,y, suffix = c(".x", ".y"), distance = FALSE) {
  hits <- hits_nearest_left(x, y)
  distance_col <- set_distance_col(distance)
  expand_by_hits(x, y, suffix, hits, hits_mcols_to_keep = distance_col)
}

#' @rdname ranges-nearest
#' @export
join_nearest_right <- function(x, y, suffix = c(".x", ".y"), distance = FALSE) {
  hits <- hits_nearest_right(x, y)
  distance_col <- set_distance_col(distance)
  expand_by_hits(x, y, suffix, hits, hits_mcols_to_keep = distance_col)
}

#' @rdname ranges-nearest
#' @export
join_nearest_upstream <- function(x, y,  suffix = c(".x", ".y"), distance = FALSE) { UseMethod("join_nearest_upstream")}

#' @export
join_nearest_upstream.GenomicRanges <- function(x, y,  suffix = c(".x", ".y"), distance = FALSE) {
  hits <- hits_nearest_upstream(x, y)
  distance_col <- set_distance_col(distance)
  expand_by_hits(x, y, suffix, hits, hits_mcols_to_keep = distance_col)
}

#' @rdname ranges-nearest
#' @importFrom IRanges nearest
#' @export
join_nearest_downstream <- function(x, y,  suffix = c(".x", ".y"), distance = FALSE) { UseMethod("join_nearest_downstream")}

#' @export
join_nearest_downstream.GenomicRanges <- function(x, y, suffix = c(".x", ".y"), distance = FALSE) {
  hits <- hits_nearest_downstream(x, y)
  distance_col <- set_distance_col(distance)
  expand_by_hits(x, y, suffix, hits, hits_mcols_to_keep = distance_col)
}

# -----------------------
# Directed hits helpers

hits_nearest <- function(x, y){ 
  UseMethod("hits_nearest")
}
hits_nearest_left <- function(x, y){ 
  UseMethod("hits_nearest_left")
}
hits_nearest_right <- function(x, y){ 
  UseMethod("hits_nearest_right")
}
hits_nearest_upstream <- function(x, y){ 
  UseMethod("hits_nearest_upstream")
}
hits_nearest_downstream <- function(x, y){ 
  UseMethod("hits_nearest_downstream")
}

#' @importFrom IRanges distanceToNearest
hits_nearest.IRanges <- function(x, y){
  make_hits(x, y, distanceToNearest, select = "arbitrary")
}

#' @importFrom GenomicRanges distanceToNearest
hits_nearest.GenomicRanges <- function(x, y){
  make_hits(x, y, distanceToNearest, select = "arbitrary", ignore.strand = TRUE)
}

#' @importFrom IRanges distanceToNearest
hits_nearest_left.IRanges <- function(x, y){
  hits <- make_hits(x, y, distanceToNearest, select = "all")
  get_hits_left(x, y, hits)
}

#' @importFrom GenomicRanges distanceToNearest
hits_nearest_left.GenomicRanges <- function(x, y){
  hits <- make_hits(x, y, distanceToNearest, select = "all", ignore.strand = TRUE)
  get_hits_left(x,y, hits)
}

#' @importFrom IRanges distanceToNearest
hits_nearest_right.IRanges <- function(x, y){
  hits <- make_hits(x, y, distanceToNearest, select = "all")
  get_hits_right(x, y, hits)
}

#' @importFrom GenomicRanges distanceToNearest
hits_nearest_right.GenomicRanges <- function(x, y){
  hits <- make_hits(x, y, distanceToNearest, select = "all", ignore.strand = TRUE)
  get_hits_right(x, y, hits)
}

#' @importFrom GenomicRanges distanceToNearest
hits_nearest_upstream.GenomicRanges <- function(x, y){
  hits <- distanceToNearest(x,y, select = "all", ignore.strand = FALSE)
  get_hits_upstream(x, y, hits)
}

#' @importFrom GenomicRanges distanceToNearest
hits_nearest_downstream.GenomicRanges <- function(x, y){
  hits <- distanceToNearest(x,y, select = "all", ignore.strand = FALSE)
  get_hits_downstream(x, y, hits)
}


# -----------------------
# Hits logic helpers

get_hits_left <- function(x, y, hits){
  mcols(hits)$is_left <- start(y[subjectHits(hits)]) <= start(x[queryHits(hits)]) &
    end(y[subjectHits(hits)]) <= start(x[queryHits(hits)])
  hits[mcols(hits)$is_left]
}

get_hits_right <- function(x, y, hits){
  mcols(hits)$is_right <- end(x[queryHits(hits)]) <= start(y[subjectHits(hits)])
  hits[mcols(hits)$is_right]
}

get_hits_upstream <- function(x, y, hits){
  mcols(hits)$is_right <- end(x[queryHits(hits)]) <= start(y[subjectHits(hits)])
  mcols(hits)$is_left <- start(x[queryHits(hits)]) >= start(y[subjectHits(hits)])
  mcols(hits)$direction <- strand(x[queryHits(hits)])
  mcols(hits)$is_upstream <- ifelse(mcols(hits)$direction == "+",
                                    mcols(hits)$is_left,
                                    mcols(hits)$is_right)
  hits[mcols(hits)$is_upstream]
}

get_hits_downstream <- function(x, y, hits){
  mcols(hits)$is_right <- end(x[queryHits(hits)]) <= start(y[subjectHits(hits)])
  mcols(hits)$is_left <- start(x[queryHits(hits)]) >= start(y[subjectHits(hits)])
  mcols(hits)$direction <- strand(x[queryHits(hits)])
  # on positive strand nearest downtream will be on the right
  # on negative strand nearest downstream will be on the left
  mcols(hits)$is_downstream <- ifelse(mcols(hits)$direction == "-",
                                    mcols(hits)$is_left,
                                    mcols(hits)$is_right)

  hits[mcols(hits)$is_downstream]
}

set_distance_col <- function(distance){
  if (distance == TRUE){
    keep_distance <- "distance"
  } else if (is.character(distance) & length(distance) == 1) {
    keep_distance <- c("distance")
    names(keep_distance) <- distance
  } else if (distance == FALSE) {
    keep_distance = NULL
  } else {
    stop("Invalid argument passed to distance")
  }
  
  return(keep_distance)
}
