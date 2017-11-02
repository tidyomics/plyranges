#' @export
tidyr::nest

#' Nest a GenomicRanges object
#'
#' @param data A Ranges object
#' @param ... A selection of metadata column variables to nest by.
#' @param .key The name of the resulting nested column.
#'
#' @importFrom tidyr nest
#' @importFrom IRanges splitAsList
#' @importFrom tidyselect vars_select
#' @importFrom GenomicRanges granges
#' @rdname ranges-nest
#' @method unnest GenomicRanges
#' @export
nest.GenomicRanges <- function(data, ..., .key = ".data") {

  dots <- quos(...)
  dot_names <- unlist(Map(function(.) quo_name(.), dots))

  if (length(dot_names) == 0L) {
    stop("Unable to fully nest a GenomicRanges object.")
  }

  selected_mcols <- tidyselect::vars_select(names(mcols(data)), !!! dots)

  groups_mcols <- setdiff(names(mcols(data)), selected_mcols)
  print(groups_mcols)
  groups_rle <- as(lapply(groups_mcols, function(.) Rle(mcols(data)[[.]])), "RleList")

  groups_split <- splitAsList(mcols(data)[, selected_mcols], groups_rle)

  nest_rng <- unique(granges(data))
  mcols(nest_rng) <- do.call("DataFrame",
                             lapply(groups_mcols, function(.) unique(mcols(data)[[.]])))
  names(mcols(nest_rng)) <- groups_mcols
  mcols(nest_rng)[[.key]] <- as(groups_split, "DataFrame")
  nest_rng

}


#' @export
tidyr::unnest

#' Expand list-columns in a Ranges object
#'
#' @param data A Ranges object
#' @param ... list-column names to unnest
#' @param .drop not implemented
#' @param .id not implemented
#' @param .sep not implemented
#'
#' @importFrom tidyr unnest
#' @method nest GenomicRanges
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

  mcols_list <- as(mcols(data)[, which_unnest], "List")
  which_unnest <- names(mcols_list)
  # this is slow
  rle_inx <-  lapply(mcols_list, lengths)
  # this is wrong if there a different groupings within each list-col
  rle_inx <- Rle(seq_along(data), Reduce(unique, rle_inx))

  # also slow
  list_vals <- Map(function(.) unlist(mcols_list[[.]]), which_unnest)
  list_vals <- do.call("DataFrame", list_vals)

  expand_rng <- granges(data)[rle_inx]

  mcols(expand_rng)  <- mcols(data)[rle_inx, !list_cols_pos]
  names(mcols(expand_rng)) <- names(list_cols_pos[!list_cols_pos])
  mcols(expand_rng) <- cbind(mcols(expand_rng), list_vals)
  expand_rng
}
