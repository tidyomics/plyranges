#' @import dplyr
#' @importFrom tidyselect all_of
mutate_mcols <- function(df, dots) {
  
  core_cols <- intersect(colnames(df),c("start", "end", "width", "seqnames", "strand"))
  
  df <- mutate(df, !!!dots) %>%
    ungroup() %>%
    dplyr::select(-all_of(core_cols)) %>%
    as("DataFrame")
}

#' @importFrom methods selectMethod
mutate_core <- function(.data, dots) {
  
  overscope <- overscope_ranges(.data)
  .mutated <- overscope_eval_update(overscope, dots)
  
  for (col in names(.mutated)) {
    modifier <- match.fun(paste0("set_", col))
    .data <- modifier(.data, .mutated[[col]])
  }
  .data
}



mutate_rng <- function(.data, ...) {
  
  dots <- set_dots_named(...)
  
  col_names <- names(dots)
  core_cols <- col_names %in% c("start", "end", "width", "seqnames", "strand")
  if (any(col_names %in% "")) {
    stop("mutate must have name-variable pairs as input", call. = FALSE)
  }
  
  if(any(core_cols)) {
    .data <- mutate_core(.data, dots[core_cols])
  }
  
  if(any(!core_cols)) {
    mcols(.data) <- mutate_mcols(as.data.frame(.data),
                          dots[!core_cols])
  }
  return(.data)
}

#' @importFrom dplyr group_by
mutate_grp <- function(.data, ...) {
  dots <- set_dots_named(...)
  
  
  col_names <- names(dots)
  core_cols <- col_names %in% c("start", "end", "width", "seqnames", "strand")
  if (any(col_names %in% "")) {
    stop("mutate must have name-variable pairs as input", call. = FALSE)
  }
  
  if(any(core_cols)) {
    inx <- .group_rows(.data)
    rng <- unname(S4Vectors::split(.data@delegate, .data@group_indices))
    rng <- S4Vectors::endoapply(rng, function(x) {
      mutate_core(x, dots[core_cols])  
    })
    
    rng <- unlist(rng)[BiocGenerics::order(unlist(inx))]
    new(class(.data),
        delegate = rng, 
        group_keys =  .data@group_keys, 
        group_indices = .data@group_indices,
        n = .data@n )
  }
  
  if (any(!core_cols)) {
    grps <- groups(.data)
    
    df <- group_by(as.data.frame(.data),!!!grps)

    mcols(.data) <-
      mutate_mcols(df,
                   dots[!core_cols])
  }
  return(.data)
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
  mutate_rng(.data, ...)
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