#' Group a Ranges by one or more variables
#'
#' @param .data a Ranges object
#' @param ... Variable names to group by. These can be either metadata columns
#' or the core variables of a Ranges.
#'
#' @description \code{group_by} takes a Ranges object and converts it
#' to a RangesGrouped object. Operations are performed on these objects by
#' their group. \code{ungroup} removes grouping.
#'
#' @importFrom dplyr group_by
#' @importFrom rlang quo_name quos syms
#' @importFrom methods new
#' @method group_by GRanges
#' @return a \code{GroupedRanges} object
#' @name group_by-ranges
#' @rdname group_by-ranges
#' @export
#' @examples
#' df <- data.frame(start = 1:10, width = 5,  seqnames = "seq1",
#' strand = sample(c("+", "-", "*"), 10, replace = TRUE), gc = runif(10))
#' rng <- Ranges(df)
#' rng
#' # grouping does not change appearance or shape or Ranges
#' rng_by_strand <- rng %>% group_by(strand)
#' rng_by_strand
#' # grouping with other verbs
#' rng_by_strand %>% summarise(gc = mean(gc))
#' rng_by_strand %>% filter(gc == min(gc))
#' # remove grouping with ungroup
#' rng_by_strand %>%
#'   ungroup() %>%
#'   summarise(gc = mean(gc))
#'
#'
group_by.GRanges <- function(.data, ...) {
  capture_groups <- quos(...)
  groups <- lapply(capture_groups, function(x) quo_name(x))
  groups <- syms(groups)
  new("GRangesGrouped", .data, groups = groups)

}


#' @param x a Ranges object
#'
#' @importFrom dplyr ungroup
#' @rdname group_by-ranges
#' @method ungroup GRangesGrouped
#' @export
ungroup.GRangesGrouped <- function(x, ...) {
  groups <- as.character(unlist(groups(x)))
  capture_groups <- quos(...)
  ungroups <- lapply(capture_groups, function(x) quo_name(x))
  if (length(ungroups) == 0L) {
    GRanges(x)
  } else {
    groups_update <- setdiff(groups, ungroups)
    groups_update <- syms(groups_update)
    new("GRangesGrouped", x, groups = groups_update)
  }
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

#' @rdname group_by-ranges
#' @method ungroup IRangesGrouped
#' @export
ungroup.IRangesGrouped <- function(x, ...) {
  groups <- as.character(unlist(groups(x)))
  capture_groups <- quos(...)
  ungroups <- lapply(capture_groups, function(x) quo_name(x))
  if (length(ungroups) == 0L) {
    IRanges(x)
  } else {
    groups_update <- setdiff(groups, ungroups)
    groups_update <- syms(groups_update)
    new("IRangesGrouped", x, groups = groups_update)
  }
}



#' Extract groupings from a RangesGrouped object
#' @param x a RangesGrouped object
#' @importFrom dplyr groups
#' @method groups GRangesGrouped
#' @rdname groups-ranges
#' @export
groups.GRangesGrouped <- function(x) { x@groups }

#' @importFrom dplyr group_vars
#' @method group_vars GRangesGrouped
#' @rdname groups-ranges
#' @export
group_vars.GRangesGrouped <- function(x) { as.character(unlist(x@groups)) }

#' @method groups GRangesGrouped
#' @rdname groups-ranges
#' @export
groups.IRangesGrouped <- function(x) { x@groups }

#' @importFrom dplyr group_vars
#' @method group_vars IRangesGrouped
#' @rdname groups-ranges
#' @export
group_vars.IRangesGrouped <- function(x) { as.character(unlist(x@groups)) }

# returns groups as split as GRangesList or RangesList
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
    groups <- as.character(unlist(groups))
    df_groups <- unique(as.data.frame(rng_list)[, groups, drop = FALSE])
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
