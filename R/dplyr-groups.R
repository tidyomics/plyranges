
#' Group a Ranges by one or more variables
#'
#' @param .data a Ranges object
#' @param ... Variable names to group by. These can be either metadata columns
#' or the core variables of a Ranges.
#'
#' @importFrom dplyr group_by
#' @importFrom rlang quo_name quos syms
#' @return a \code{GroupedGRanges} object
#' @export
#' @rdname group_by
group_by.GRanges <- function(.data, ...) {
  capture_groups <- quos(...)
  groups <- lapply(capture_groups, function(x) quo_name(x))
  groups <- syms(groups)
  new("GRangesGrouped", .data, groups = groups)

}

#' @rdname group_by
group_by.IRanges <- function(.data, ...) {
  capture_groups <- quos(...)
  groups <- lapply(capture_groups, function(x) quo_name(x))
  groups <- syms(groups)
  new("IRangesGrouped", .data, groups = groups)
}

#' @importFrom dplyr groups
groups.GRangesGrouped <- function(x) { x@groups }
groups.IRangesGrouped <- function(x) { x@groups }

# returns groups as split as GRangesList or RangesList
split_groups <- function(.data_grouped) {
  groups <- groups(.data_grouped)
  rng_env <- as.env(.data_grouped, parent.frame())
  rle_groups <- lapply(groups, eval_bare, env = rng_env)
  rle_groups <- as(rle_groups, "RleList")
  S4Vectors::splitAsList(.data_grouped, rle_groups)
}
