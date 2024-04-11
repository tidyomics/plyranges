mutate_mcols <- function(.data, dots) {
  
  overscope <- overscope_ranges(.data)
  .mutated <- overscope_eval_update(overscope, dots)
  
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


#' @import dplyr
#' @importFrom tidyselect all_of
mutate_mcols_grp <- function(df, dots) {
  
  core_cols <- intersect(colnames(df),c("start", "end", "width", "seqnames", "strand"))
  
  df <- mutate(df, !!!dots) %>%
    ungroup() %>%
    dplyr::select(-tidyselect::all_of(core_cols)) %>%
    as("DataFrame")
}

#' @importFrom methods selectMethod
mutate_core <- function(.data, dots) {
  
  overscope <- overscope_ranges(.data)
  .mutated <- overscope_eval_update(overscope, dots)
  
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



mutate_rng <- function(.data, ...) {
  
  dots <- set_dots_named(...)
  
  col_names <- names(dots)
  core_cols <- col_names %in% c("start", "end", "width", "seqnames", "strand")
  if (any(col_names %in% "")) {
    stop("mutate must have name-variable pairs as input", call. = FALSE)
  }
  
  if(any(core_cols)) {
    .data <- mutate_core(.data, dots)
  }
  
  if(any(!core_cols)) {
    .data <- mutate_mcols(.data,
                          dots)
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
      mutate_core(x, dots)  
    })
    
    rng <- unlist(rng)[BiocGenerics::order(unlist(inx))]
    new(class(.data),
        delegate = rng, 
        group_keys =  .data@group_keys, 
        group_indices = .data@group_indices,
        n = .data@n )
  }
  
  if (any(!core_cols)) {
    grps <- group_vars(.data)

    df <- as.data.frame(ungroup(.data))

    df <- group_by(df,!!!rlang::syms(grps))

    mcols(.data) <-
      mutate_mcols_grp(df,
                   dots)
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