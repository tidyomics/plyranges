summarize_rng <- function(.data, ...) {
  dots <- UQS(...)
  gr_env <- as.env(.data, parent.frame())
  overscope <- new_overscope(gr_env, parent.env(gr_env))
  on.exit(overscope_clean(overscope))
  summarised <- vector("list", length(dots))
  names(summarised) <- names(dots)

  for (i in seq_along(summarised)) {
    summarised[[i]] <- overscope_eval_next(overscope, dots[[i]])
    # sometimes we want to compute on previously constructed columns
    # we can do this by binding the evaluated expression to
    # the overscope environment
    new_col <- names(dots)[[i]]
    rlang::env_bind(overscope$.env, !!new_col := summarised[[i]])
  }
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
  split_ranges <- split_groups(.data, populate_mcols = TRUE, drop = TRUE)
  groups_summary <- lapply(split_ranges, summarize_rng, dots)
  groups_summary <- lapply(groups_summary, as_tibble)
  groups_df <- as_tibble(as.data.frame(mcols(split_ranges)))
  groups_summary <- bind_rows(groups_summary)
  bind_cols(groups_df, groups_summary)

}
