summarize_rng <- function(.data, dots) {
  overscope <- overscope_ranges(.data)
  ans <- overscope_eval_update(overscope, dots, FALSE)
  
  # maintain list columns instead of collapsing them
  is_list <- vapply(ans, function(.) is(., "List") || is(., "list"), logical(1))
  
  # compress list columns
  if (any(is_list)) {
    nr <- check_n(.data)
    for (i in which(is_list)) {
      if (length(ans[[i]]) == 1) {
        ans[[i]] <- as(rep(ans[[i]], nr), "CompressedList")
      } else {
        
        if (all(lengths(ans[[i]]) == length(ans[[i]][[1]]))) {
          ans[[i]] <- as(BiocGenerics::Reduce(S4Vectors::pc, ans[[i]]), 
                         "CompressedList")
          # check length is equal to number of rows or records
          stopifnot(length(ans[[i]]) == nr)
        }
        

      }
    } 
  }
  
  results <- DataFrame(ans)
  rownames(results) <- NULL
  results
}

check_n <- function(.data) {
  if (is(.data, "GroupedGenomicRanges") || is(.data, "GroupedIntegerRanges")) {
    return(.data@n)
  }
  1L
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
  dots <- set_dots_named(...)
  summarize_rng(.data, dots)
}

#' @method summarise DelegatingGenomicRanges
#' @export
summarise.DelegatingGenomicRanges <- function(.data, ...) {
  dots <- set_dots_named(...)
  delegate <- .data@delegate
  summarize_rng(delegate, dots)
}

#' @method summarise DelegatingGenomicRanges
#' @export
summarise.DelegatingIntegerRanges <- function(.data, ...) {
  dots <- set_dots_named(...)
  delegate <- .data@delegate
  summarize_rng(delegate, dots)
}

#' @importFrom rlang !!! enquos
#' @importFrom dplyr bind_cols bind_rows
#' @method summarise GroupedGenomicRanges
#' @export
summarise.GroupedGenomicRanges <- function(.data, ...) {
  dots <- set_dots_named(...)
  md <- .data@group_keys
  res <- cbind(md, summarize_rng(.data, dots))
  res[order(res[, group_vars(.data)]), ]
}

#' @method summarise GroupedIntegerRanges
#' @export
summarise.GroupedIntegerRanges <- summarise.GroupedGenomicRanges