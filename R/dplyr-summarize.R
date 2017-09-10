summarize_rng <- function(.data, ...) {
  dots <- UQS(...)

  gr_env <- as.env(.data, parent.frame())
  overscope <- new_overscope(gr_env, parent.env(gr_env))
  on.exit(overscope_clean(overscope))
  summarised <- lapply(dots, overscope_eval_next, overscope = overscope)
  summarised
}

#' Aggregate a GenomicRanges
#'
#' @param .data a GRanges object
#' @param ... Name-value pairs of summary functions.
#'
#' @return a \link[tibble]{tibble}
#' @importFrom tibble as_tibble
#' @importFrom dplyr summarize
#' @importFrom dplyr summarise
summarise.GRanges <- function(.data, ...) {

  dots <- quos(...)
  as_tibble(summarize_rng(.data, dots))

}


#' @importFrom rlang UQS quos
#' @importFrom dplyr bind_cols bind_rows
summarise.GRangesGrouped <- function(.data, ...) {

  dots <- quos(...)

  split_ranges <- split_groups(.data, populate_mcols = TRUE)
  groups_summary <- bind_rows(lapply(split_ranges, summarize_rng, dots))
  groups_df <- as.data.frame(mcols(split_ranges))
  as_tibble(bind_cols(groups_df, groups_summary))

}

group_levels <- function(x) {
  if (is.factor(unique(x))) {
    levels(x)
  } else {
    unique(x)
  }
}
