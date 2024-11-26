# ported from HelloRanges
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

zero_cols <- function(x) {
  empty <- new(class(S4Vectors::decode(x)))
  if (is(x, "Rle")) return(Rle(values = empty))
  empty
}

na_ranges <- function(x) {
  na <- na_cols(x)
  if (length(mcols(x)) > 0) {
    mcols(na) <- DataFrame(lapply(mcols(x), na_cols))
  }
  na
}

na_dframe <- function(dframe, nrows) {
  
  empty <- new("DFrame", nrows = as.integer(nrows))
  
  # IRanges special case for mcols is null
  if (is.null(ncol(dframe))) return(NULL)
  if (ncol(dframe) == 0) return(empty)
  
  tform <- na_cols
  if (nrows == 0) tform <- zero_cols
  
  for (i in colnames(dframe)) {
      empty[[i]] <- tform(dframe[[i]])
  }

  empty
  
}


add_na_seqlevels <- function(x) {
  na <- na_cols(granges(x))
  seqlevels(x) <- union(seqlevels(na), seqlevels(x))
  x
}

.join_overlap_left <- function(x, y, suffix, f_in, distance = FALSE, ...) {
  # generate hits
  hits <- make_hits(x, y, f_in, ...)
  # overlaps found
  left <- x[queryHits(hits), ]
  right <- y[subjectHits(hits), ]
  mcols(left) <- mcols_overlaps_update(left, right, suffix)

  # overlaps not found
  only_left <- rep(TRUE, queryLength(hits))
  only_left[queryHits(hits)] <- FALSE

  # ranges object
  rng_only_left <- x[only_left]
  
  # if there are no right mcols create an empty DataFrame
  # otherwise propagate missing values to the outer frame
  mcols_outer <- na_dframe(mcols(right), sum(only_left))
  
  if (!is.null(mcols(rng_only_left))) {
    mcols(rng_only_left) <- cbind(mcols(rng_only_left), mcols_outer)

  } else {
    mcols(rng_only_left) <- mcols_outer
  }

  names(mcols(rng_only_left)) <- names(mcols(left))

  left_outer <- c(left, rng_only_left)
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
join_overlap_left.IntegerRanges <- function(x,y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  .join_overlap_left(x,y, suffix, findOverlaps, 
                 maxgap = maxgap, 
                 minoverlap = minoverlap, 
                 type = "any",
                 select = "all")
}

#' @export
join_overlap_left.GenomicRanges <- function(x, y, maxgap = -1L, minoverlap = 0L,
                                            suffix = c(".x", ".y"), distance = FALSE) {
  .join_overlap_left(x,
                     y, 
                     suffix, 
                     findOverlaps, 
                     maxgap = maxgap, 
                     minoverlap = minoverlap, 
                     type = "any",
                     select = "all",
                     ignore.strand = TRUE,
                     distance = distance)
}

#' @rdname overlap-joins
#' @export
join_overlap_left_within <- function(x, y, maxgap, minoverlap, suffix = c(".x", ".y")) {
  UseMethod("join_overlap_left_within")
}

#' @export
join_overlap_left_within.IntegerRanges <- function(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  .join_overlap_left(x,y, suffix, findOverlaps, 
                     maxgap = maxgap, 
                     minoverlap = minoverlap, 
                     type = "within",
                     select = "all")
}

#' @export
join_overlap_left_within.GenomicRanges <- function(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  .join_overlap_left((x),
                     (y), 
                     suffix, 
                     findOverlaps, 
                     maxgap = maxgap, 
                     minoverlap = minoverlap, 
                     type = "within",
                     select = "all",
                     ignore.strand = TRUE)
}

#' @rdname overlap-joins
#' @export
join_overlap_left_directed <- function(x, y, maxgap, minoverlap, suffix = c(".x", ".y")) {
  UseMethod("join_overlap_left_directed")
}


#' @export
join_overlap_left_directed.GenomicRanges <- function(x, y, maxgap = -1L, minoverlap = 0L,
                                                     suffix = c(".x", ".y"), distance = FALSE) {
  .join_overlap_left((x),
                     (y), 
                     suffix, 
                     findOverlaps, 
                     maxgap = maxgap, 
                     minoverlap = minoverlap, 
                     type = "any",
                     select = "all",
                     ignore.strand = FALSE,
                     distance = distance)
}

#' @rdname overlap-joins
#' @export
join_overlap_left_within_directed <- function(x, y, maxgap, minoverlap, suffix = c(".x", ".y")) {
  UseMethod("join_overlap_left_within_directed")
}


#' @export
join_overlap_left_within_directed.GenomicRanges <- function(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y")) {
  .join_overlap_left(x,
                     y, 
                     suffix, 
                     findOverlaps, 
                     maxgap = maxgap, 
                     minoverlap = minoverlap, 
                     type = "within",
                     select = "all",
                     ignore.strand = FALSE)
}
