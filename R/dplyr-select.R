select_rng <- function(.data, ...) {
  dots <- UQS(...)
  col_names <- unlist(lapply(dots, quo_name))
  col_syms <- syms(col_names)
  names(col_syms) <- col_names

  rng_os <- overscope_ranges(.data)
  on.exit(overscope_clean(rng_os))

  rng_df <- lapply(col_syms, overscope_eval_next, overscope = rng_os)
  rng_df <- as_tibble(as.data.frame(rng_df))

  rng <- try(Ranges(rng_df), silent = TRUE)

  if (is(rng, "try-error")) {
    return(rng_df)
  } else {
    return(rng)
  }
}


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
#' @rdname ranges-select
#' @method select GenomicRanges
#' @export
select.GenomicRanges <- function(.data, ...) {

  dots <- quos(...)

  # no selection? - return .data
  if (length(dots) == 0) {
    return(.data)
  }

  # if the quo is named return an error
  stopifnot(length(names(dots)) != 0)

  select_rng(.data, dots)

}

#' @rdname ranges-select
#' @method select Ranges
#' @export
select.Ranges <- function(.data, ...) {
  dots <- quos(...)

  # no selection? - return .data
  if (length(dots) == 0) {
    return(.data)
  }

  # if the quo is named return an error
  stopifnot(length(names(dots)) != 0)

  select_rng(.data, dots)
}
