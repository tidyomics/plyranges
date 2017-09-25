#' Select parts of the GRanges by name
#'
#' @param .data a \code{Ranges} object
#' @param ... One or more column names including the core components
#' of the Ranges object - if possible select will preserve the Ranges object,
#' however in cases where this is not possible select will return a tibble.
#' @details Note that if a core component of a Ranges is dropped or selected
#' without the other required components (this includes the seqnames, strand, start, end,
#' width names), then select will return a tibble.
#' @note select doesn't currently support slicing syntax or integer based selection.
#' @return a Ranges object or a tibble
#' @seealso \link[dplyr]{select}
#' @importFrom dplyr select
select.GenomicRanges <- function(.data, ...) {

  dots <- quos(...)

  # no selection? - return .data
  if (length(dots) == 0) {
    return(.data)
  }

  # if the quo is named
  stopifnot(length(names(dots)) != 0)

  col_names <- unlist(lapply(dots, quo_name))
  col_syms <- syms(col_names)
  names(col_syms) <- col_names

  rng_os <- overscope_ranges(.data)

  rng_df <- lapply(col_syms, overscope_eval_next, overscope = rng_os)
  rng_df <- as_tibble(as.data.frame(rng_df))

  rng <- try(Ranges(rng_df), silent = TRUE)

  if (is(rng, "try-error")) {
    return(rng_df)
  } else {
    return(rng)
  }

  # col_slice <- grep(":", col_names)
  #
  # if (length(col_slice) > 0) {
  #   start_col <- sub("^(.*?):(.*)",
  #                    "\\1",
  #                    col_names[col_slice])
  #   end_col <- sub("^(.*?):(.*)", "\\2", col_names[col_slice])
  #   col_slice_index <- lapply(seq_along(col_slice),
  #                             function(i) {
  #                               start <- grep(start_col[i], mcol_names)
  #                               end <- grep(end_col[i], mcol_names)
  #                               seq.int(start, end, by = 1)
  #                               })
  #   col_slice_index <- unlist(col_slice_index)
  #   col_names <- col_names[-col_slice]
  # } else {
  #   col_slice_index <- integer(0)
  # }
  #
  # if (length(col_names) > 0) {
  #   print(col_names)
  #   col_index <- vapply(col_names,
  #                       function(x)  {
  #                         grep(sub("^-", "", x), mcol_names) *
  #                           ifelse(grepl("^-", x), -1L, 1L)
  #                       },
  #                       integer(1))
  #   col_drops <- sign(col_index) == -1
  # } else {
  #   col_index <- integer(0)
  # }
  #
  # if (any(col_drops)) {
  #   .data <- .data[, col_index[col_drops]]
  #   col_index <- col_index[!col_drops]
  # }
  #
  # .data[, c(col_slice_index,col_index)]

}
