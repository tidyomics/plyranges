#' @method select DeferredGenomicRanges
#' @importFrom Rsamtools bamWhat<- bamTag<-
#' @export
select.DeferredGenomicRanges <- function(.data, ..., .drop_ranges = FALSE) {
 
  # check to see if we need to update cache, by selecting columns
  # that are not present in the delegate
  vars <- tbl_vars(.data@delegate)
  names(vars) <- vars 
  pos <- try(tidyselect::eval_select(rlang::expr(c(...)),
                                      vars,
                                      exclude = parallelVectorNames(.data@delegate)),
              silent = TRUE)
  
  
  if (is_empty_delegate(.data) || is(pos, "try-error")) {
    .data@ops <- select(.data@ops, ..., .drop_ranges = .drop_ranges)
    .data@delegate <- load_genomic_file(.data@ops)
  } else {
    .data@delegate <- select(.data@delegate, ..., .drop_ranges = .drop_ranges)
  }
  .data
}

#' @method filter DeferredGenomicRanges
#' @importFrom Rsamtools scanBamFlag bamFlag<-
#' @export
filter.DeferredGenomicRanges <- function(.data, ...) {
  dots <- quos(...)
  # check to see if we need to override cache 
  # this happens if any valid_flag_filters are available
  dot_names <- gsub("^!", "", unlist(lapply(dots, quo_name)))
  update_cache <- any(dot_names %in% unlist(names(valid_flag_filters())))
  # clear the cache
  if (update_cache & !is_empty_delegate(.data)) {
    .data@delegate <- GRanges()
  }
  
  if (is_empty_delegate(.data)) {
    .data@ops <- filter(.data@ops, ...)
    .data@delegate <- load_genomic_file(.data@ops)
  } else {
    .data@delegate <- filter(.data@delegate, ...)
  }
  .data
}



#' @importFrom Rsamtools bamWhich<-
#' @export
filter_by_overlaps.DeferredGenomicRanges <- function(x, y, maxgap = -1L, minoverlap = 0L) {
  if (is_empty_delegate(x)) {
    x@ops <- filter_by_overlaps(x@ops, y, maxgap, minoverlap)
    x@delegate <- load_genomic_file(x@ops)
  } else {
    x@delegate <- filter_by_overlaps(x@delegate, y, maxgap, minoverlap)
  }
  x
}

#' @method summarise DeferredGenomicRanges
#' @export
summarise.DeferredGenomicRanges <- function(.data, ...) {
  summarize(load_delegate(.data), ...)
}

#' @method mutate DeferredGenomicRanges
#' @export
mutate.DeferredGenomicRanges <- function(.data, ...) {
  mutate(load_delegate(.data), ...)
}

#' @method group_by DeferredGenomicRanges
#' @export
group_by.DeferredGenomicRanges <- function(.data, ...) {
  group_by(load_delegate(.data), ...)
}
