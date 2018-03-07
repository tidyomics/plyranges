#' Interweave a pair of Ranges objects together
#'
#' @param left,right Ranges objects.
#' @param .id When supplied a new column that represents the origin column
#' and is linked to each row of the resulting Ranges object.
#'
#' @details The output of `interweave()` takes pairs of Ranges
#' objects and combines them into a single  Ranges object. If an .id
#' argument is supplied, an origin column with name .id is created indicated which side
#' the resulting Range comes from (eit)
#'
#' @examples
#'gr <- as_granges(data.frame(start = 10:15,
#'                             width = 5,
#'                             seqnames = "seq1",
#'                             strand = c("+", "+", "-", "-", "+", "*")))
#' interweave(flank_left(gr, width = 5L), flank_right(gr, width = 5L))
#' interweave(flank_left(gr, width = 5L), flank_right(gr, width = 5L), .id = "origin")
#' @importFrom S4Vectors zipup Pairs
#' @export
#' @rdname ranges-interweave
interweave <- function(left, right, .id = NULL) {
  if (!is.null(.id)) {
    stopifnot(is.character(.id) & length(.id) == 1L)
    mcols(left)[[.id]] <- "left"
    mcols(right)[[.id]] <- "right"
  }
  unlist(zipup(Pairs(left, right)))
}


#' Combine Ranges by concatentating them together
#'
#' @param ... Ranges objects to combine. Each argument can be a Ranges object,
#' or a list of Ranges objects.
#' @param .id Ranges object identifier. When .id is supplied a new column
#' is created that links each row to the original Range object. The contents
#' of the column correspond to the named arguments or the names of the list
#' supplied.
#'
#' @importFrom rlang dots_values have_name
#' @rdname ranges-bind
#' @examples
#'gr <- as_granges(data.frame(start = 10:15,
#'                             width = 5,
#'                             seqnames = "seq1"))
#'gr2 <- as_granges(data.frame(start = 11:14,
#'                             width = 1:4,
#'                             seqnames = "seq2"))
#'bind_ranges(gr, gr2)
#'
#'bind_ranges(a = gr, b = gr2, .id = "origin")
#'
#'bind_ranges(gr, list(gr, gr2), gr2)
#'
#'bind_ranges(list(a = gr, b = gr2), c = gr, .id = "origin")
#'
#' @export
bind_ranges <- function(..., .id = NULL) {
  x <- unlist(rlang::dots_values(...))

  if (all(vapply(x, function(x) is(x, "GenomicRanges"), logical(1)))) {
    to_class <- "GRangesList"
  } else if (all(vapply(x, function(x) is(x, "Ranges"), logical(1)))) {
    to_class <- "IRangesList"
  } else {
    stop("Cannot bind objects of different classes together, ... must be
         all Ranges or GenomicRanges object.", call. = FALSE)
  }

  if (!is.null(.id)) {
    stopifnot(is.character(.id) & length(.id) == 1L)
    any_empty <-  vapply(x, function(x) length(x) == 0, logical(1))
    any_noname <- rlang::have_name(x)
    if (!all(any_empty | any_noname)) {
      x <- Filter(length, x)
      names(x) <- seq_along(x)
    }
  }


  if (!is.null(.id)) {
    rng <- unlist(as(x, to_class), use.names = TRUE)
    mcols(rng)[[.id]] <-  names(rng)
    names(rng) <- NULL
    return(rng)
  } else {
    return(unlist(as(x, to_class), use.names = FALSE))
  }
}
