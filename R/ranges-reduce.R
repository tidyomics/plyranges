# ranges-reduce
reduce_rng <- function(.data, reduced, dots) {

  os <- S4Vectors::as.env(.data, 
                          parent.frame(), 
                          tform = function(col) {
                            unname(IRanges::extractList(col, mcols(reduced)$revmap))
                          })
  os <- new_data_mask(os, top = parent.env(os))
  

  mcols(reduced) <- DataFrame(overscope_eval_update(os, dots, FALSE))
  return(reduced)
}

reduce_single <- function(.data, ..., rfun = reduce) {
  dots <- set_dots_named(...)
  if (length(dots) == 0L) {
    return(rfun(.data))
  }
  reduced <- rfun(.data, with.revmap = TRUE)
  reduce_rng(.data, reduced, dots)
}

reduce_by_grp <- function(.data, ..., rfun = IRanges::reduce) {
  dots <- set_dots_named(...)
  delegate <- .data@delegate
  inx <- .data@inx
  groups <- group_vars(.data)
  split_ranges <- extractList(delegate, inx)
  mcols(split_ranges) <- mcols(inx)
  
  if (length(dots) == 0L) {
    rng <- IRanges::stack(rfun(split_ranges))
    if (is(rng, "IntegerRanges")) {
      mcols(rng) <- mcols(inx[mcols(rng)[["name"]]])[, groups, drop = FALSE]
    } else {
      mcols(rng) <- mcols(rng)[, groups, drop = FALSE]
    }
    return(rng)
  }
  
  reduced <- rfun(split_ranges, with.revmap = TRUE)
  rng <- Map(function(x, y) {reduce_rng(x, y, dots)}, 
             split_ranges, 
             reduced)
  rng <- as(rng, "List")
  mcols(rng) <- mcols(split_ranges)
  rng <- IRanges::stack(rng)
  # drop 'name' col from stack
  mcols(rng) <- mcols(rng)[, -1, drop = FALSE]
  # reorder columns
  mcols(rng) <- mcols(rng)[, rev(seq_along(mcols(rng)))]
  rng
}

#' Reduce then aggregate a Ranges object
#'
#' @param .data a Ranges object to reduce
#' @param ... Name-value pairs of summary functions.
#'
#' @return a Ranges object with the
#' @rdname ranges-reduce
#' @importFrom IRanges reduce
#' @importFrom utils relist
#' @examples
#' set.seed(10)
#' df <- data.frame(start = sample(1:10), 
#'                  width = 5,  
#'                  seqnames = "seq1",
#'                  strand = sample(c("+", "-", "*"), 10, replace = TRUE), 
#'                  gc = runif(10))
#'                  
#' rng <- as_granges(df)
#' rng %>% reduce_ranges()
#' rng %>% reduce_ranges(gc = mean(gc))
#' rng %>% reduce_ranges_directed(gc = mean(gc))
#' 
#' x <- data.frame(start = c(11:13, 2, 7:6), 
#'                width=3, 
#'                id=sample(letters[1:3], 6, replace = TRUE),
#'                score= sample(1:6))
#' x <- as_iranges(x)
#' x %>% reduce_ranges()
#' x %>% reduce_ranges(score = sum(score))
#' x %>% group_by(id) %>% reduce_ranges(score = sum(score))
#' @export
reduce_ranges <- function(.data, ...) { UseMethod("reduce_ranges") }

#' @method reduce_ranges IntegerRanges
#' @export
reduce_ranges.IntegerRanges <- function(.data, ...) {
  reduce_single(.data, ...)
}

#' @method reduce_ranges GroupedIntegerRanges
#' @export
reduce_ranges.GroupedIntegerRanges <- function(.data, ...) {
  reduce_by_grp(.data, ...)
}


#' @method reduce_ranges GroupedGenomicRanges
#' @export
reduce_ranges.GroupedGenomicRanges <- function(.data, ...) {
  reduce_by_grp(.data, ..., 
                rfun = function(x, ...) {
                  reduce(x, ..., ignore.strand = TRUE)
                  })
}

#' @method reduce_ranges GenomicRanges
#' @export
reduce_ranges.GenomicRanges <- function(.data, ...) {
  reduce_single(.data, ..., 
                rfun = function(x, ...) {
                  reduce(x, ..., ignore.strand = TRUE)
                })
}

#' @rdname ranges-reduce
#' @export
reduce_ranges_directed <- function(.data, ...) {
  UseMethod("reduce_ranges_directed")
}

#' @importFrom IRanges reduce
#' @method reduce_ranges_directed GenomicRanges
#' @export
reduce_ranges_directed.GenomicRanges <- function(.data, ...) {
  reduce_single(.data, ..., 
                rfun = function(x, ...) {
                  reduce(x, ..., ignore.strand = FALSE)
                })
}

#' @method reduce_ranges_directed GroupedGenomicRanges
#' @export
reduce_ranges_directed.GroupedGenomicRanges <- function(.data, ...) {
  reduce_by_grp(.data, ..., 
                rfun = function(x, ...) {
                  reduce(x, ..., ignore.strand = FALSE)
                })
}