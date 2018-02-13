# ported from HelloRanges
na_ranges <- function(x) {
  na <- na_cols(x)
  if (length(mcols(x)) > 0) {
    mcols(na) <- DataFrame(lapply(mcols(x), na_cols))
  }
  na
}

na_cols <- function(x) {
  if (is(x, "IRanges"))
    return(IRanges(0L, -1L))
  if (is(x, "GRanges"))
    return(GRanges(".", IRanges(0L, -1L)))
  if (is(x, "Rle"))
    x <- S4Vectors::decode(x)
  ans <- x[NA_integer_]
  if (is(x, "Rle"))
    ans <- S4Vectors::Rle(ans)
  ans
}

add_na_seqlevels <- function(x) {
  na <- na_cols(granges(x))
  seqlevels(x) <- union(seqlevels(na), seqlevels(x))
  x
}

outer_left_join <- function(x, y, hits, suffix) {
  # overlaps found
  left <- x[queryHits(hits), ]
  right <- y[subjectHits(hits), ]
  mcols(left) <- mcols_overlaps_update(left, right, suffix)

  # overlaps not found
  only_left <- rep(TRUE, queryLength(hits))
  only_left[queryHits(hits)] <- FALSE

  # ranges object
  rng_only_left <- x[only_left]
  mcols_outer <- rep(DataFrame(lapply(mcols(right), na_cols)),
                              sum(only_left))
  if (!is.null(mcols(rng_only_left))) {
    mcols(rng_only_left) <- cbind(mcols(rng_only_left), mcols_outer)

  } else {
    mcols(rng_only_left) <- mcols_outer
  }

  names(mcols(rng_only_left)) <- names(mcols(left))


  left_outer <-   c(left, rng_only_left)
  left_outer <- left_outer[order(c(queryHits(hits), which(only_left)))]
  left_outer
}

#' @importFrom S4Vectors queryLength decode
#' @importFrom GenomeInfoDb seqlevels seqlevels<-
#' @rdname overlap-joins
#' @export
join_overlap_left <- function(x, y, maxgap, minoverlap, suffix = c(".x", ".y")) {
  UseMethod("join_overlap_left")
}

#' @export
join_overlap_left.Ranges <- function(x,y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  hits <- findOverlaps(x,y, maxgap, minoverlap, type = "any", select = "all")
  outer_left_join(x,y,hits, suffix)
}

#' @export
join_overlap_left.GenomicRanges <- function(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  x <- add_na_seqlevels(x)
  y <- add_na_seqlevels(y)
  hits <- findOverlaps(x,y, maxgap, minoverlap,
                       type = "any", select = "all", ignore.strand = TRUE)

  outer_left_join(x,y,hits, suffix)
}

#' @rdname overlap-joins
#' @export
join_overlap_left_within <- function(x, y, maxgap, minoverlap, suffix = c(".x", ".y")) {
  UseMethod("join_overlap_left_within")
}

#' @export
join_overlap_left_within.Ranges <- function(x, y, maxgap, minoverlap, suffix = c(".x", ".y")) {
  hits <- findOverlaps(x,y, maxgap, minoverlap, type = "within", select = "all")
  outer_left_join(x, y, hits, suffix)
}

#' @export
join_overlap_left_within.GenomicRanges <- function(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  x <- add_na_seqlevels(x)
  y <- add_na_seqlevels(y)

  hits <- findOverlaps(x,y, maxgap, minoverlap,
                       type = "within", select = "all", ignore.strand = TRUE)

  outer_left_join(x,y,hits, suffix)
}

#' @rdname overlap-joins
#' @export
join_overlap_left_directed <- function(x, y, maxgap, minoverlap, suffix = c(".x", ".y")) {
  UseMethod("join_overlap_left_directed")
}


#' @export
join_overlap_left_directed.GenomicRanges <- function(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  x <- add_na_seqlevels(x)
  y <- add_na_seqlevels(y)

  hits <- findOverlaps(x,y, maxgap, minoverlap,
                       type = "any", select = "all", ignore.strand = FALSE)
  outer_left_join(x,y,hits,suffix)
}
