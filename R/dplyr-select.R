select_rng <- function(.data, .drop_ranges, ...) {
  
  var_names <- tbl_vars(.data)
  names(var_names) <- var_names
  if (.drop_ranges) {
    pos <- tidyselect::eval_select(rlang::expr(c(...)), var_names)
    
    .env <- overscope_ranges(.data)
    
    ans <- lapply(syms(names(pos)), eval_tidy, data = .env)
    names(ans) <- names(pos)
    return(as(ans, "DataFrame"))
  } else {
    
    core <- S4Vectors::parallelVectorNames(.data)
    var_names <- var_names[!(var_names %in% core)]

    
    pos <- try(tidyselect::eval_select(rlang::expr(c(...)), var_names),
               silent = TRUE)
    
    if (is(pos, "try-error") || !length(var_names)) {
      invalid_cols <- paste(core, collapse = ", ")
      stop(paste0("Cannot select/rename the following columns: ",invalid_cols))
    }
  }
  
  mcols(.data) <- mcols(.data)[ , pos, drop = FALSE]
  names(mcols(.data)) <- names(pos)
  
  .data
}


#' Select metadata columns of the Ranges object by name or position
#'
#' @param .data a `Ranges` object
#' @param ... One or more metadata column names.
#' @param .drop_ranges If TRUE select will always return a tibble. In this
#' case, you may select columns that form the core part of the Ranges object.
#' @details Note that by default select only acts on the metadata columns (and
#' will therefore return a Ranges object) if a core component of a Ranges is dropped or selected
#' without the other required components (this includes the seqnames, strand, start, end,
#' width names), then select will throw an error unless .drop_ranges is set to TRUE.
#' @return a Ranges object or a tibble
#' @seealso [dplyr::select()]
#' @importFrom dplyr select
#' @importFrom tidyselect eval_select
#' @importFrom S4Vectors parallelVectorNames
#' @examples
#' df <- data.frame(start = 1:10, width = 5,  seqnames = "seq1",
#' strand = sample(c("+", "-", "*"), 10, replace = TRUE), gc = runif(10), counts = rpois(10, 2))
#' rng <- as_granges(df)
#' select(rng, -gc)
#' select(rng, gc)
#' select(rng, counts, gc)
#' select(rng, 2:1)
#' select(rng, seqnames, strand, .drop_ranges = TRUE)
#' @rdname ranges-select
#' @method select Ranges
#' @export
select.Ranges <- function(.data, ..., .drop_ranges = FALSE) {
  dots <- rlang::enquos(...)
  # no selection? - return .data
  if (length(dots) == 0) {
    return(.data)
  }

  select_rng(.data, .drop_ranges, ...)
}

#' @method select DelegatingGenomicRanges
#' @export
select.DelegatingGenomicRanges <- function(.data, ..., .drop_ranges = FALSE) {
  dots <- rlang::enquos(...)
  # no selection? - return .data
  if (length(dots) == 0) {
    return(.data)
  }

  delegate <- .data@delegate
  delegate <- select_rng(delegate, .drop_ranges,  dots)
  if (.drop_ranges) {
    return(delegate)
  } else {
    .data@delegate <- delegate
  }
  .data
}

#' @method select DelegatingIntegerRanges
#' @export
select.DelegatingIntegerRanges <- select.DelegatingGenomicRanges