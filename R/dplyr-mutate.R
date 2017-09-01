mutate_mcols <- function(.data, .mutated) {
  all_cols <- names(.mutated)
  only_mcols <- !(all_cols %in%
                    c("start", "end", "width", "seqnames", "strand"))
  .mutated <- .mutated[only_mcols]
  update_cols <- all_cols[only_mcols]

  matches_mcols <- match(update_cols, names(mcols(.data)))
  idx_mcols <- !is.na(matches_mcols)

  if (any(idx_mcols)) {
    mcols(.data)[matches_mcols[idx_mcols]] <- .mutated[idx_mcols]
  }

  if (!all(idx_mcols)) {
    mcols(.data) <- do.call("DataFrame",
                            list(mcols(.data), .mutated[!idx_mcols]))
  }

  .data
}

mutate_core <- function(.data, .mutated) {
  all_cols <- names(.mutated)
  core_cols <- all_cols[all_cols %in%
                           c("start", "end", "width", "seqnames", "strand")]
  if (length(core_cols == 0)) {
    .data
  }

  for (col in core_cols) {
    accessor <- selectMethod(paste0(col, "<-"), class(.data))
    .data <- accessor(.data, value = .mutated[[col]])
  }

  .data
}

#' modify a GenomicRanges object
#'
#' @param .data a \code{GRanges} object
#' @param ... dots!
#'
#' @importFrom dplyr mutate
#' @return a GRanges object
mutate.GRanges <- function(.data, ...) {
  dots <- UQS(quos(...))

  print(dots)
  col_names <- names(dots)
  if (any(col_names %in% "")) {
     stop("mutate must have name-variable pairs as input", .call = FALSE)
  }

  gr_env <- as.env(.data, parent.frame())
  overscope <- new_overscope(gr_env, parent.env(gr_env))

  .mutated <- lapply(dots, overscope_eval_next, overscope = overscope)
  on.exit(overscope_clean(overscope))

  .data <- mutate_core(.data, .mutated)
  mutate_mcols(.data, .mutated)

}

mutate.IRanges <- function(.data, ...) {

}
