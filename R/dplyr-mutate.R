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
      mcols(.data) <- S4Vectors::DataFrame(.mutated[!idx_mcols])
    } else {
      mcols(.data) <- S4Vectors::DataFrame(list(mcols(.data), 
                                                .mutated[!idx_mcols]))
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
    modifier <- match.fun(paste0("set_", col))
    .data <- modifier(.data, .mutated[[col]])
  }
  .data
}


mutate_rng <- function(.data, dots) {
  col_names <- names(dots)
  if (any(col_names %in% "")) {
    stop("mutate must have name-variable pairs as input", call. = FALSE)
  }

  overscope <- overscope_ranges(.data)
  .mutated <- overscope_eval_update(overscope, dots)
  .data <- mutate_core(.data, .mutated)
  mutate_mcols(.data, .mutated)
}

mutate_grp <- function(.data, ...) {
  dots <- set_dots_named(...)
  delegate <- .data@delegate
  inx <- .data@inx
  rng <- lapply(inx, function(i) {
    x <- delegate[i]
    mutate_rng(x, dots)    
  })
  rng <- bind_ranges(rng)[order(unlist(inx))]
  names(rng) <- NULL
  new(class(.data),
      delegate = rng, 
      groups =  groups(.data), 
      inx = inx)
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
#' # mutate can be used in conjuction with anchoring to resize ranges
#' rng %>%
#'     mutate(width = 10)
#' # by default width modfication fixes by start
#' rng %>%
#'     anchor_start() %>%
#'     mutate(width = 10)
#' # fix by end or midpoint
#' rng %>%
#'     anchor_end() %>%
#'     mutate(width = width + 1)
#' rng %>%
#'     anchor_center() %>%
#'     mutate(width = width + 1)
#' # anchoring by strand
#' rng %>%
#'     anchor_3p() %>%
#'     mutate(width = width * 2)
#' rng %>%
#'     anchor_5p() %>%
#'     mutate(width = width * 2)
#' @export
mutate.Ranges <- function(.data, ...) {
  dots <- set_dots_named(...)
  mutate_rng(.data, dots)
}

#' @method mutate AnchoredIntegerRanges
#' @export
mutate.AnchoredIntegerRanges <- mutate.Ranges

#' @method mutate AnchoredGenomicRanges
#' @export
mutate.AnchoredGenomicRanges <- mutate.Ranges

#' @method mutate DelegatingGenomicRanges
#' @export
mutate.DelegatingGenomicRanges <- function(.data, ...) {
  dots <- set_dots_named(...)
  delegate <- .data@delegate
  .data@delegate <- mutate_rng(delegate, dots)
  return(.data)
}

#' @method mutate DelegatingIntegerRanges
#' @export
mutate.DelegatingIntegerRanges <- mutate.DelegatingGenomicRanges

#' @method mutate GroupedGenomicRanges
#' @export
mutate.GroupedGenomicRanges <- function(.data, ...) {
  mutate_grp(.data, ...)
}

#' @method mutate GroupedIntegerRanges
#' @export
mutate.GroupedIntegerRanges <- mutate.GroupedGenomicRanges