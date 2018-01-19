summarize_rng <- function(.data, ...) {
  dots <- UQS(...)
  overscope <- overscope_ranges(.data)
  on.exit(overscope_clean(overscope))
  overscope_eval_update(overscope, dots)
}


#' Aggregate a Ranges object
#'
#' @param .data a Ranges object
#' @param ... Name-value pairs of summary functions.
#'
#' @return a \link[S4Vectors]{DataFrame}
#' @importFrom S4Vectors rbind cbind
#' @importFrom dplyr summarise summarize
#' @method summarise GenomicRanges
#' @rdname ranges-summarise
#' @examples
#' df <- data.frame(start = 1:10, width = 5,  seqnames = "seq1",
#' strand = sample(c("+", "-", "*"), 10, replace = TRUE), gc = runif(10))
#' rng <- as_granges(df)
#' rng %>% summarise(gc = mean(gc))
#' rng %>% group_by(strand) %>% summarise(gc = mean(gc))
#' @export
summarise.GenomicRanges <- function(.data, ...) {

  dots <- quos(...)

  DataFrame(summarize_rng(.data, dots))

}

#' @method summarise Ranges
#' @rdname ranges-summarise
#' @export
summarise.Ranges <- function(.data, ...) {

  dots <- quos(...)
  DataFrame(summarize_rng(.data, dots))

}

#' @importFrom rlang UQS quos
#' @importFrom dplyr bind_cols bind_rows
#' @method summarise GRangesGrouped
#' @rdname ranges-summarise
#' @export
summarise.GRangesGrouped <- function(.data, ...) {

  dots <- quos(...)
  split_ranges <- split_groups(.data, populate_mcols = TRUE, drop = TRUE)
  groups_summary <- lapply(split_ranges, summarize_rng, dots)
  groups_summary <- do.call(rbind, lapply(groups_summary, as, "DataFrame"))
  cbind(mcols(split_ranges), groups_summary)

}

#' @method summarise IRangesGrouped
#' @rdname ranges-summarise
#' @export
summarise.IRangesGrouped <- function(.data, ...) {
  dots <- quos(...)
  split_ranges <- split_groups(.data, populate_mcols = TRUE, drop = TRUE)
  groups_summary <- lapply(split_ranges, summarize_rng, dots)
  groups_summary <- do.call(rbind, lapply(groups_summary, as, "DataFrame"))
  cbind(mcols(split_ranges), groups_summary)
}

