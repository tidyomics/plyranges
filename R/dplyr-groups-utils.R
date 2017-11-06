#' Split a Grouped Ranges object into a RangesList
#' @param .data_grouped A grouped ranges object
#' @param populate_mcols add grouping to metadata slot of RangesList
#' @param drop drop unused group levels
#'
#' @importFrom methods as
#' @importFrom S4Vectors Rle
#' @importFrom rlang eval_bare
split_groups <- function(.data_grouped, populate_mcols = FALSE, drop = TRUE) {
  groups <- groups(.data_grouped)
  rng_env <- as.env(.data_grouped, parent.frame())
  list_groups <- lapply(groups, function(x) {
    grp <- eval_bare(x, env = rng_env)
    as(grp, "Rle")
  })

  names(list_groups) <- as.character(groups)

  rle_groups <- as(list_groups, "RleList")
  rng_list <- IRanges::splitAsList(.data_grouped, rle_groups, drop = drop)
  if (populate_mcols) {
    group_vars <- group_vars(.data_grouped)
    if (is(rng_list, "GRangesList")) {
      df_groups <- unique(as.data.frame(rng_list)[, group_vars, drop = FALSE])
    } else {
      df_groups <- lapply(rng_list, function(x) {
        slice <- cbind(as.data.frame(x), mcols(x))
        unique(slice[, group_vars, drop = FALSE])
      })
      df_groups <- do.call("rbind", df_groups)
    }
    rownames(df_groups) <- NULL
    mcols(rng_list) <- df_groups
  }
  rng_list
}

group_levels <- function(x) {
  if (is.factor(x)) {
    unique(as.character(x))
  } else {
    unique(x)
  }
}
