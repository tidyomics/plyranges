summarize_rng <- function(.data, dots) {
  overscope <- overscope_ranges(.data)
  on.exit(overscope_clean(overscope))
  DataFrame(overscope_eval_update(overscope, dots))
}

rename_dots <- function(dots) {
  unnamed_dots <- nchar(names(dots)) == 0
  rename_dots <- unlist(lapply(dots[unnamed_dots], rlang::quo_name))
  names(dots)[unnamed_dots] <- rename_dots
  dots
}

#' Aggregate a Ranges object
#'
#' @param .data a Ranges object
#' @param ... Name-value pairs of summary functions.
#'
#' @return a [S4Vectors::DataFrame()]
#' @seealso [dplyr::summarise()]
#' @importFrom S4Vectors rbind cbind
#' @importFrom dplyr summarise summarize
#' @rdname ranges-summarise
#' @examples
#' df <- data.frame(start = 1:10, width = 5,  seqnames = "seq1",
#' strand = sample(c("+", "-", "*"), 10, replace = TRUE), gc = runif(10))
#' rng <- as_granges(df)
#' rng %>% summarise(gc = mean(gc))
#' rng %>% group_by(strand) %>% summarise(gc = mean(gc))
#' @method summarise Ranges
#' @rdname ranges-summarise
#' @export
summarise.Ranges <- function(.data, ...) {
  dots <- quos(...)
  dots <- rename_dots(dots)
  summarize_rng(.data, dots)
}

#' @method summarise DelegatingGenomicRanges
#' @rdname ranges-summarise
#' @export
summarise.DelegatingGenomicRanges <- function(.data, ...) {
  dots <- quos(...)
  dots <- rename_dots(dots)
  delegate <- .data@delegate
  summarize_rng(delegate, dots)
}

#' @method summarise DelegatingGenomicRanges
#' @rdname ranges-summarise
#' @export
summarise.DelegatingIntegerRanges <- function(.data, ...) {
  dots <- quos(...)
  dots <- rename_dots(dots)
  delegate <- .data@delegate
  summarize_rng(delegate, dots)
}

#' @importFrom rlang UQS quos
#' @importFrom dplyr bind_cols bind_rows
#' @method summarise GroupedGenomicRanges
#' @rdname ranges-summarise
#' @export
summarise.GroupedGenomicRanges <- function(.data, ...) {
  dots <- quos(...)
  dots <- rename_dots(dots)
  delegate <- .data@delegate
  inx <- .data@inx
  groups_summary <- as(
    lapply(
          extractList(delegate, inx),
          summarize_rng, 
          dots
          ),
          "List"
  )
  cbind(mcols(inx), unlist(groups_summary, use.names = FALSE))
}

#' @method summarise GroupedIntegerRanges
#' @rdname ranges-summarise
#' @export
summarise.GroupedIntegerRanges <- summarise.GroupedGenomicRanges