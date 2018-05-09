#' @method select DeferredGenomicRanges
#' @importFrom Rsamtools bamWhat<- bamTag<-
#' @export
select.DeferredGenomicRanges <- function(.data, ..., .drop_ranges = FALSE) {
  if (is_empty_delegate(.data)) {
    .data@ops <- select(.data@ops, ..., .drop_ranges = .drop_ranges)
  } else {
    .data@delegate <- select(.data@delegate, ..., .drop_ranges = .drop_ranges)
  }
  .data
}

#' @method filter DeferredGenomicRanges
#' @importFrom Rsamtools scanBamFlag bamFlag<-
#' @export
filter.DeferredGenomicRanges <- function(.data, ...) {
  if (is_empty_delegate(.data)) {
    .data@ops <- filter(.data@ops, ...)
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
