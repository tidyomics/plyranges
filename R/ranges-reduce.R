# ranges-reduce
group_by_revmap <- function(.data, indexes, groups) {
  group_vars <- syms(c(group_vars(.data), "revmap"))
  .data <- ungroup(.data)
  .data <- .data[unlist(indexes)]
  mcols(.data)[["revmap"]] <- groups
  return(group_by(.data, !!!group_vars))
} 

make_key_rle <- function(x) {
  Rle(as.integer(S4Vectors::runValue(x)), S4Vectors::runLength(x))
}

make_revmap_rle <- function(x) {
  Rle(seq_along(x), elementNROWS(x))
}

unsplice <- function(x, y) x[unlist(y)]



add_revmap_grouping <- function(.data, key, revmap) {
  # lookup indexes for each group
  lookup <- S4Vectors::splitAsList(revmap, key)
  inx <- .group_rows(.data)
  # indexes we will expand by
  indexes <- S4Vectors::mendoapply(unsplice, inx, lookup)
  # groupings from revmap
  groups <- make_revmap_rle(revmap)
  group_by_revmap(.data, indexes, groups)
}

reduce_single <- function(.data, ..., rfun = reduce) {
  dots <- set_dots_named(...)
  if (length(dots) == 0L) {
    return(rfun(.data))
  }

  reduced <- rfun(.data, with.revmap = TRUE)
  
  .data <- group_by_revmap(.data,
                           mcols(reduced)[["revmap"]],
                           make_revmap_rle(mcols(reduced)[["revmap"]]))
  
  sd <- summarise(.data, !!!dots)
  
  sd <- sd[order(sd[["revmap"]]), -which(names(sd) == "revmap")[1], drop = FALSE]
  
  mcols(reduced) <- sd
  reduced
}

reduce_by_grp <- function(.data, ..., rfun = IRanges::reduce) {
  dots <- set_dots_named(...)

  by_groups <- dplyr::group_split(.data)
  
  if (length(dots) == 0L) {
    rng <- IRanges::stack(rfun(by_groups))
    sd <- dplyr::group_keys(.data)
    key <- make_key_rle(mcols(rng)[["name"]])
    mcols(rng) <- sd[key, , drop = FALSE]
    return(rng)
  }
  
  rng <- IRanges::stack(rfun(by_groups, with.revmap = TRUE))

  .data <- add_revmap_grouping(.data, 
                               mcols(rng)[["name"]], 
                               mcols(rng)[["revmap"]])
  
  sd <- summarise(.data, !!!dots)
  sd <- sd[order(sd[["revmap"]]), -which(names(sd) == "revmap"), drop = FALSE]
  
  mcols(rng) <- sd
  rng
}

#' Reduce then aggregate a Ranges object
#'
#' @param .data a Ranges object to reduce
#' @param min.gapwidth Ranges separated by a gap of at least min.gapwidth positions are not merged. 
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
reduce_ranges <- function(.data, min.gapwidth = 1L, ...) { UseMethod("reduce_ranges") }

#' @method reduce_ranges IntegerRanges
#' @export
reduce_ranges.IntegerRanges <- function(.data, min.gapwidth = 1L, ...) {
  reduce_single(.data, ..., 
                rfun = function(x, ...) {
                  reduce(x, ..., min.gapwidth = min.gapwidth)
                })
}

#' @method reduce_ranges GroupedIntegerRanges
#' @export
reduce_ranges.GroupedIntegerRanges <- function(.data, min.gapwidth = 1L, ...) {
  reduce_by_grp(.data, ..., 
                rfun = function(x, ...) {
                  reduce(x, ..., min.gapwidth = min.gapwidth)
                })
}

#' @method reduce_ranges GroupedGenomicRanges
#' @export
reduce_ranges.GroupedGenomicRanges <- function(.data, min.gapwidth = 1L, ...) {
  reduce_by_grp(.data, ..., 
                rfun = function(x, ...) {
                  reduce(x, ..., min.gapwidth = min.gapwidth, ignore.strand = TRUE)
                  })
}

#' @method reduce_ranges GenomicRanges
#' @export
reduce_ranges.GenomicRanges <- function(.data, min.gapwidth = 1L, ...) {
  reduce_single(.data, ..., 
                rfun = function(x, ...) {
                  reduce(x, ..., min.gapwidth = min.gapwidth, ignore.strand = TRUE)
                })
}

#' @rdname ranges-reduce
#' @export
reduce_ranges_directed <- function(.data, min.gapwidth = 1L, ...) {
  UseMethod("reduce_ranges_directed")
}

#' @importFrom IRanges reduce
#' @method reduce_ranges_directed GenomicRanges
#' @export
reduce_ranges_directed.GenomicRanges <- function(.data, min.gapwidth = 1L, ...) {
  reduce_single(.data, ..., 
                rfun = function(x, ...) {
                  reduce(x, ..., min.gapwidth = min.gapwidth, ignore.strand = FALSE)
                })
}

#' @method reduce_ranges_directed GroupedGenomicRanges
#' @export
reduce_ranges_directed.GroupedGenomicRanges <- function(.data, min.gapwidth = 1L, ...) {
  reduce_by_grp(.data, ..., 
                rfun = function(x, ...) {
                  reduce(x, ..., min.gapwidth = min.gapwidth, ignore.strand = FALSE)
                })
}
