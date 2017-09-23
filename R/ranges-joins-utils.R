# ranges-joins-utils.R
ranges_vars <- function(x) {
  x_env <- as.env(x, parent.frame())
  vars_rng <-ls(x_env)
  vars_rng <- vars_rng[!(vars_rng %in% "names")]
  vars_mcols <- ls(parent.env(x_env))
  c(vars_rng, vars_mcols)
}

# dplyr's join syntax uses a function called tbl_vars to get
# variable names, using this function will enable a Ranges to be copied through
# as a data.frame in a join.
tbl_vars.GenomicRanges <- function(x) {
  ranges_vars(x)
}


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
