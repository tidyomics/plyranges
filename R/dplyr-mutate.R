#' modify a GenomicRanges object
#'
#' @param .data a \code{GRanges} object
#' @param ... dots!
#'
#' @return a GRanges object
mutate.GRanges <- function(.data, ...) {
  dots <- quos(...)

  col_names <- names(dots)

  invalid_cols <- any(grepl("strand|start|end|width", col_names))
  if (invalid_cols) {
    stop("Invalid selection, only metadata columns may be mutated.")
  }
  gr_env <- as.env(.data, parent.frame())
  overscope <- new_overscope(gr_env, parent.env(gr_env))

  transformed <- lapply(dots, overscope_eval_next, overscope = overscope)
  on.exit(overscope_clean(overscope))

  update_cols <- names(transformed)

  # mcols mutate
  matches_mcols <- match(update_cols, names(mcols(.data)))
  idx_mcols <- !is.na(matches_mcols)

  if (any(idx_mcols)) {
    mcols(.data)[matches_mcols[idx_mcols]] <- transformed[idx_mcols]
  }

  if (!all(idx_mcols)) {
    mcols(.data) <- do.call("DataFrame",
                            list(mcols(.data), transformed[!idx_mcols]))
  }

  .data
}
