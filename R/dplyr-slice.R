generate_slice <- function(.data, dots) {
  os <- overscope_ranges(.data)
  inx <- overscope_eval_update(os, dots, bind_envir = FALSE)
  check_num <- vapply(inx, is.numeric, logical(1))
  
  if (any(!check_num)) {
    stop("slice condition does not evaluate to an integer or numeric vector")
  }
  inx
}

#' Choose rows by their position
#' 
#' @param .data a `Ranges` object
#' @param ... Integer row values indicating rows to keep. If `.data` is
#' a [`GroupedGenomicRanges`] or [`GroupedIntegerRanges`] object.
#' @param .preserve when FALSE (the default) the grouping structure is 
#' recomputed, otherwise it is kept as is. Currently ignored. 
#' @return a GRanges object
#' @method slice Ranges
#' @importFrom dplyr slice
#' @rdname slice-ranges
#'
#' @examples
#' df <- data.frame(start = 1:10,
#'                  width = 5,
#'                  seqnames = "seq1",
#'                  strand = sample(c("+", "-", "*"), 10, replace = TRUE),
#'                  gc = runif(10))
#' rng <- as_granges(df)
#' dplyr::slice(rng, 1:2)
#' dplyr::slice(rng, -n())
#' dplyr::slice(rng, -5:-n())
#' 
#' by_strand <- group_by(rng, strand)
#' 
#' # slice with group by finds positions within each group
#' slice(by_strand, n())
#' slice(by_strand, which.max(gc))
#' 
#' # if the index is beyond the number of groups slice are ignored
#' slice(by_strand, 1:3)
#' 
#' 
#' @export
slice.Ranges <- function(.data, ..., .preserve = FALSE) {
  dots <- set_dots_unnamed(...)
  inx <- generate_slice(.data, dots)
  inx <- Reduce(union, inx)
  .data[inx, ]
}

#' @rdname slice-ranges
#' @method slice GroupedGenomicRanges
#' @export
slice.GroupedGenomicRanges <- function(.data, ..., .preserve = FALSE) {
  dots <- set_dots_unnamed(...)
  inx <- .group_rows(.data)
  rng <- .data@delegate
  mch <- as(
    lapply(inx, 
           function(i) unlist(generate_slice(rng[i], dots))
    ),
    "IntegerList"
  )
  inx <- inx[mch]
  inx <- inx[!is.na(inx)]
  
  rng <- rng[sort(unlist(inx))]
  
  group_by(rng, !!!groups(.data))
  
}

#' @rdname slice-ranges
#' @method slice GroupedIntegerRanges
#' @export
slice.GroupedIntegerRanges <- slice.GroupedGenomicRanges