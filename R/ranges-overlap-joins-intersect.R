.join_intersect <- function(x, y, suffix, f_in = findOverlapPairs, ...) {

  pairs <- f_in(x, y, ...)
  strand_arg <- quos(...)
  strand_arg <- strand_arg[names(strand_arg) == "ignore.strand"]
  if (length(strand_arg) == 1L) {
      inner_rng <- pintersect(pairs, ignore.strand = eval_tidy(strand_arg[[1]]))
  } else {
    inner_rng <- pintersect(pairs)
  }
  left_rng <- first(pairs)
  right_rng <- second(pairs)
  mcols(inner_rng) <- mcols_overlaps_update(left_rng, right_rng, suffix)
  inner_rng
}

#' Join by overlapping Ranges
#'
#' @param x,y Objects representing ranges
#' @param maxgap,minoverlap The maximimum gap between intervals as an integer
#' greater than or equal to zero. The minimum amount of overlap between intervals
#' as an integer greater than zero, accounting for the maximum gap.
#' @param suffix Character to vectors to append to common columns in x and y
#' (default = `c(".x", ".y")`).
#'
#' @details The function [join_overlap_intersect()] finds
#' the genomic intervals that are the overlapping ranges between x and y and
#' returns a new ranges object with metadata columns from x and y.
#'
#' The function [join_overlap_inner()] is equivalent to [find_overlaps()].
#'
#' The function [join_overlap_left()] performs a left outer join between x
#' and y. It returns all ranges in x that overlap or do not overlap ranges in y
#' plus metadata columns common to both. If there is no overlapping range
#' the metadata column will contain a missing value.
#'
#' The function [join_overlap_self()] find all overlaps between a ranges
#' object x and itself.
#'
#' All of these functions have two suffixes that modify their behavior.
#' The `within` suffix, returns only ranges in x that are completely
#' overlapped within in y. The `directed` suffix accounts for the strandedness
#' of the ranges when performing overlaps.
#'
#' @return a GRanges object
#'
#' @importFrom S4Vectors first second DataFrame
#' @importFrom IRanges findOverlapPairs
#' @examples
#' x <- as_iranges(data.frame(start = c(11, 101), end = c(21, 201)))
#' y <- as_iranges(data.frame(start = c(10, 20, 50, 100, 1),
#'                            end = c(19, 21, 105, 202, 5)))
#' 
#' # self
#' join_overlap_self(y)
#' 
#' # intersect takes common interval
#' join_overlap_intersect(x,y)
#'
#' # within
#' join_overlap_intersect_within(x,y)
#' 
#' # left, and inner join, it's often useful having an id column here
#' y <- y %>% mutate(id = 1:n()) 
#' x <- x %>% mutate(id = 1:n()) 
#' join_overlap_inner(x,y)
#' join_overlap_left(y,x, suffix = c(".left", ".right"))
#' 
#' @rdname overlap-joins
#' @export
#' @seealso [join_overlap_self()], [join_overlap_left()], [find_overlaps()]
join_overlap_intersect <- function(x, y, maxgap, minoverlap, suffix = c(".x", ".y")) {
  UseMethod("join_overlap_intersect")
}

#' @export
join_overlap_intersect.IntegerRanges <- function(x, y, maxgap = -1L, minoverlap = 0L,
                                          suffix = c(".x", ".y")) {
  .join_intersect(x, y, 
                  suffix, 
                  type = "any", 
                  maxgap = maxgap, 
                  minoverlap = minoverlap)
}

#' @export
join_overlap_intersect.GenomicRanges <- function(x, y, maxgap = -1L, minoverlap = 0L,
                                                 suffix = c(".x", ".y")) {
  .join_intersect(x, y, 
                  suffix, 
                  type = "any", 
                  maxgap = maxgap, 
                  minoverlap = minoverlap,
                  ignore.strand = TRUE)

}
#' @rdname overlap-joins
#' @export
join_overlap_intersect_within <- function(x, y, maxgap, minoverlap,
                                      suffix = c(".x", ".y")) {
  UseMethod("join_overlap_intersect_within")
}

#' @export
join_overlap_intersect_within.IntegerRanges <- function(x, y, maxgap = -1L, minoverlap = 0L,
                                             suffix = c(".x", ".y")) {
  .join_intersect(x, y, 
                  suffix, 
                  type = "within", 
                  maxgap = maxgap, 
                  minoverlap = minoverlap)
}

#' @export
join_overlap_intersect_within.GenomicRanges <- function(x, y, maxgap = -1L, minoverlap = 0L,
                                                    suffix = c(".x", ".y")) {
  .join_intersect(x, y, 
                  suffix, 
                  type = "within", 
                  maxgap = maxgap, 
                  minoverlap = minoverlap,
                  ignore.strand = TRUE)
}


#' @rdname overlap-joins
#' @export
join_overlap_intersect_directed <- function(x, y, maxgap, minoverlap,
                                            suffix = c(".x", ".y")) {
  UseMethod("join_overlap_intersect_directed")
}

#' @export
join_overlap_intersect_directed.GenomicRanges <- function(x, y,
                                                          maxgap = -1L,
                                                          minoverlap = 0L,
                                                          suffix = c(".x", ".y")) {
  .join_intersect(x, y, suffix,
                  type = "any",
                  maxgap = maxgap,
                  minoverlap = minoverlap,
                  ignore.strand = FALSE)
}

#' @rdname overlap-joins
#' @export
join_overlap_intersect_within_directed <- function(x, y, maxgap, minoverlap,
                                            suffix = c(".x", ".y")) {
  UseMethod("join_overlap_intersect_within_directed")
}

#' @export
join_overlap_intersect_within_directed.GenomicRanges <- function(x, y,
                                                          maxgap = -1L,
                                                          minoverlap = 0L,
                                                          suffix = c(".x", ".y")) {
  .join_intersect(x, y, suffix,
                  type = "within",
                  maxgap = maxgap,
                  minoverlap = minoverlap,
                  ignore.strand = FALSE)
}


#' @rdname overlap-joins
#' @export
join_overlap_inner <- function(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  UseMethod("find_overlaps")
}

#' @rdname overlap-joins
#' @export
join_overlap_inner_within <- function(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  UseMethod("find_overlaps_within")
}

#' @rdname overlap-joins
#' @export
join_overlap_inner_directed <- function(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  UseMethod("find_overlaps_directed")
}

#' @rdname overlap-joins
#' @export
join_overlap_inner_within_directed <- function(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  UseMethod("find_overlaps_within_directed")
}