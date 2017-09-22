# ranges-joins-utils.R
common_mcols <- function(x, y, by = NULL) {
  x_names <- names(mcols(x))
  y_names <- names(mcols(y))
  if (is.null(by)) {
    common_mcols <- intersect(x_names,  y_names)
    if (length(common_mcols) == 0) {
      stop("No common columns between x & y", call. = FALSE)
    }
    return(common_mcols)
  } else {
    named_by <- names(by)
    if (length(named_by) > 0) {
      stopifnot(named_by %in% x_names || by %in% y_names)
      by

    } else {
      stopifnot(by %in% x_names || by %in% y_names)
      by
    }
  }
}

map_common <- function(x, y, x_name, y_name) {
  if (is(x, "GRanges") & is(y, "DataFrame")) {
    intersect(mcols(x)[[x_name]], y[[y_name]])
  } else {
    intersect(mcols(x)[[x_name]], mcols(y)[[y_name]])
  }
}
index_common <- function(x, y, x_name, y_name) {
  mcols(x)[[x_name]] %in% map_common(x, y, x_name, y_name)
}

rows_common <- function(by, x, y) {
  lapply(seq_along(by),
         function(i) {
           x_name <- names(by)[[i]]
           y_name <- by[[i]]
           index_common(x, y, x_name, y_name)
         })
}

negate_rows_common <- function(by, x, y) {
  lapply(seq_along(by),
         function(i) {
           x_name <- names(by)[[i]]
           y_name <- by[[i]]
           !index_common(x, y, x_name, y_name)
         })
}

filter_common <- function(by, x, y) {
  Reduce(all, rows_common(by, x, y))
}
