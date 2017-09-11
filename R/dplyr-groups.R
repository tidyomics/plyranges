#' Group a Ranges by one or more variables
#'
#' @param .data a Ranges object
#' @param ... Variable names to group by. These can be either metadata columns
#' or the core variables of a Ranges.
#'
#' @importFrom dplyr group_by
#' @importFrom rlang quo_name quos syms
#' @importFrom methods new
#' @method group_by GRanges
#' @return a \code{GroupedRanges} object
#' @name group_by-ranges
#' @rdname group_by-ranges
#' @export
group_by.GRanges <- function(.data, ...) {
  capture_groups <- quos(...)
  groups <- lapply(capture_groups, function(x) quo_name(x))
  groups <- syms(groups)
  new("GRangesGrouped", .data, groups = groups)

}

#' @rdname group_by-ranges
#' @method group_by IRanges
#' @export
group_by.IRanges <- function(.data, ...) {
  capture_groups <- quos(...)
  groups <- lapply(capture_groups, function(x) quo_name(x))
  groups <- syms(groups)
  new("IRangesGrouped", .data, groups = groups)
}


#' Extract groupings from a RangesGrouped object
#' @param x a RangesGrouped object
#' @importFrom dplyr groups
#' @method groups GRangesGrouped
#' @rdname groups-ranges
#' @export
groups.GRangesGrouped <- function(x) { x@groups }

#' @method groups GRangesGrouped
#' @rdname groups-ranges
#' @export
groups.IRangesGrouped <- function(x) { x@groups }

# returns groups as split as GRangesList or RangesList
split_groups <- function(.data_grouped, populate_mcols = FALSE) {
  groups <- groups(.data_grouped)
  rng_env <- as.env(.data_grouped, parent.frame())
  list_groups <- lapply(groups, eval_bare, env = rng_env)
  names(list_groups) <- as.character(groups)

  rle_groups <- as(list_groups, "RleList")
  rng_list <- IRanges::splitAsList(.data_grouped, rle_groups)
  if (populate_mcols) {
    df_groups <- expand.grid(lapply(list_groups, group_levels))
    mcols(rng_list) <- df_groups
  }
  rng_list
}
