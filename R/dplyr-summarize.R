#' Aggregate a GenomicRanges
#'
#' @param .data a GRanges object
#' @param ... Name-value pairs of summary functions.
#'
#' @return a \link[tibble]{tibble}
summarise.GRanges <- function(.data, ...) {

  dots <- quos(...)

  gr_env <- as.env(.data, parent.frame())
  overscope <- new_overscope(gr_env, parent.env(gr_env))
  on.exit(overscope_clean(overscope))
  summarised <- lapply(dots, overscope_eval_next, overscope = overscope)

  as_tibble(summarised)

}


summarise.GroupedGRanges <- function(.data, ...) {

  groups <- groups(.data)
  gr_env <- as.env(.data, parent.frame())
  value <- lapply(groups, eval_bare, env = gr_env)
  value <- as(value, "RleList")
  split_ranges <- GRangesList(splitAsList(GRanges(.data), value, drop = FALSE))


  groups_data <- expand.grid(lapply(value, group_levels))
  names(groups_data) <- unlist(lapply(groups, as.character))

  by_group <- lapply(split_ranges, summarise.GRanges, !!!quos(...))
  group_by(bind_cols(groups_data, bind_rows(by_group)), !!!groups)
}

group_levels <- function(x) {
  if (is.factor(unique(x))) {
    levels(x)
  } else {
    unique(x)
  }
}
