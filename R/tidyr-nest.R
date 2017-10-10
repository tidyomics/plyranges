
#' @export
tidyr::nest

#' @importFrom tidyr nest
#' @export
nest.GenomicRanges <- function(data, ..., .key) {

}


#' @export
tidyr::unnest

#' Expand list-columns in a Ranges object
#'
#' @param data
#' @param ... list-column names to unnest
#' @param .drop not implemented
#' @param .id not implemented
#' @param .sep not implemented
#'
#' @importFrom tidyr unnest
#' @export
#' @rdname ranges-unnest
unnest.GenomicRanges <- function(data, ..., .drop, .id, .sep) {

  dots <- quos(...)
  dot_names <- unlist(Map(function(.) quo_name(.), dots))

  list_cols_pos <- unlist(Map(function(.) is(., "List"), mcols(data)))

  list_cols <- Filter(isTRUE, list_cols_pos)

  if (length(list_cols) == 0L) {
    stop("No columns to unnest.", call. = FALSE)
  }

  if (length(dot_names) == 0L) {
    which_unnest <- names(list_cols)
  } else {
    which_unnest <- intersect(dot_names, names(list_cols))
  }
  if (length(which_unnest) == 0L) {
    stop(paste("Input column(s):",
               paste0(dot_names, collapse = ","), "not found"), call. = FALSE)
  }

  # this is slow
  rle_inx <-  lapply(which_unnest, function(x) lengths(mcols(data)[[x]]))
  # this is wrong if there a different groupings within each list-col
  rle_inx <- Rle(seq_along(data), Reduce(unique, rle_inx))

  # also slow
  list_vals <- Map(function(.) unlist(mcols(data)[[.]]), which_unnest)
  list_vals <- do.call("DataFrame", list_vals)

  expand_rng <- granges(data)[rle_inx]

  mcols(expand_rng)  <- mcols(data)[rle_inx, !list_cols_pos]
  names(mcols(expand_rng)) <- names(list_cols_pos[!list_cols_pos])
  mcols(expand_rng) <- cbind(mcols(expand_rng), list_vals)
  mcols(expand_rng) <- mcols(expand_rng)[, names(list_cols_pos)]
  expand_rng
}
