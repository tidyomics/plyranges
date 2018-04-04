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
    if (is.null(mcols(.data))) {
      mcols(.data) <- do.call("DataFrame", .mutated[!idx_mcols])
    } else {
      mcols(.data) <- do.call("DataFrame",
                              list(mcols(.data), .mutated[!idx_mcols]))
    }

  }

  .data
}

#' @importFrom methods selectMethod
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


mutate_rng <- function(.data, dots) {
  #dots <- UQS(dots)
  col_names <- names(dots)
  if (any(col_names %in% "")) {
    stop("mutate must have name-variable pairs as input", .call = FALSE)
  }

  overscope <- overscope_ranges(.data)
  on.exit(overscope_clean(overscope))
  .mutated <- overscope_eval_update(overscope, dots)

  .data <- mutate_core(.data, .mutated)
  mutate_mcols(.data, .mutated)

}
#' Modify a Ranges object
#'
#' @param .data a `Ranges` object
#' @param ... Pairs of name-value expressions. The name-value pairs can either
#' create new metadata columns or modify existing ones.
#'
#' @importFrom dplyr mutate
#' @rdname mutate-ranges
#' @return a Ranges object
#' @method mutate Ranges
#'
#' @examples
#' df <- data.frame(start = 1:10,
#'                  width = 5,
#'                  seqnames = "seq1",
#'                  strand = sample(c("+", "-", "*"), 10, replace = TRUE),
#'                  gc = runif(10))
#' rng <- as_granges(df)
#'
#' # mutate adds new columns
#' rng %>%
#'     mutate(avg_gc = mean(gc), row_id = 1:n())
#' # can also compute on newly created columns
#' rng %>%
#'     mutate(score = gc * width, score2 = score + 1)
#' # group by partitions the data and computes within each group
#' rng %>%
#'     group_by(strand) %>%
#'     mutate(avg_gc = mean(gc), row_id = 1:n())
#'
#' @export
mutate.Ranges <- function(.data, ...) {
  dots <- quos(...)
  mutate_rng(.data, dots)
}

#' @rdname mutate-ranges
#' @method mutate DelegatingGenomicRanges
#' @export
mutate.DelegatingGenomicRanges <- function(.data, ...) {
  dots <- quos(...)
  delegate <- .data@delegate
  .data@delegate <- mutate_rng(delegate, dots)
  return(.data)
}

#' @rdname mutate-ranges
#' @method mutate DelegatingIntegerRanges
#' @export
mutate.DelegatingIntegerRanges <- function(.data, ...) {
  dots <- quos(...)
  delegate <- .data@delegate
  .data@delegate <- mutate_rng(delegate, dots)
  return(.data)
}

#' @rdname mutate-ranges
#' @method mutate GroupedGenomicRanges
#' @export
mutate.GroupedGenomicRanges <- function(.data, ...) {
  dots <- quos(...)
  delegate <- .data@delegate
  inx <- .data@inx
  rng <- lapply(inx, function(i) {
    x <- delegate[i]
    mutate_rng(x, dots)
    
  })

  rng <- bind_ranges(rng)[order(unlist(inx))]

  new("GroupedGenomicRanges",
      elementMetadata =  .data@elementMetadata, 
      delegate = rng, 
      groups =  groups(.data), 
      inx = inx)
}

#' @rdname mutate-ranges
#' @method mutate GroupedIntegerRanges
#' @export
mutate.GroupedIntegerRanges <- function(.data, ...) {
  dots <- quos(...)
  delegate <- .data@delegate
  inx <- .data@inx
  rng <- lapply(inx, function(i) {
    x <- delegate[i]
    mutate_rng(x, dots)  
  })

  rng <- bind_ranges(rng)[order(unlist(inx))]

  new("GroupedIntegerRanges",
      delegate = rng, 
      groups =  groups(.data), 
      inx = inx)
}
