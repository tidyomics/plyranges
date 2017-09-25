# ranges-joins-utils.R

# common columns between two Ranges
common_cols <- function(x, y, by = NULL) {

  x_names <- ranges_vars(x)
  y_names <- ranges_vars(y)

  if (is.null(by)) {
    common_cols <- intersect(x_names, y_names)
    if (length(common_mcols) == 0) {
      stop("No common columns between x & y", call. = FALSE)
    }
    return(common_cols)
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
