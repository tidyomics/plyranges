#' Group a Ranges by one or more variables
#'
#' @param .data a Ranges object.
#' @param ... Variable names to group by. These can be either metadata columns
#' or the core variables of a Ranges.
#' @param x a GroupedRanges object.
#'
#' @description The function `group_by` takes a Ranges object and defines
#' groups by one or more variables. Operations are then performed on the Ranges
#' by their "group". `ungroup()` removes grouping.
#'
#' @details
#' `group_by()` creates a new object of class `GRangesGrouped` if
#' the input is a `GRanges` object or an object of class `GroupedIntegerRanges`
#' if the input is a `IRanges` object. Both of these classes contain a slot
#' called `groups` corresponding to the names of grouping variables. They
#' also inherit from their parent classes, `Ranges` and `GenomicRanges`
#' respectively. `ungroup()` removes the grouping and will return
#' either a `GRanges` or `IRanges` object.
#'
#' @section Accessors:
#' To return grouping variables on a grouped Ranges use either
#' \itemize{
#'   \item{`groups(x)`}{Returns a list of symbols}
#'   \item{`group_vars(x)`}{Returns a character vector}
#' }
#'
#' @return The `group_by()` function will return a GroupedRanges object.
#' These have the same appearance as a regular Ranges object but with an
#' additional groups slot.
#'
#'
#' @importFrom dplyr group_by
#' @importFrom rlang quo_name quos syms
#' @importFrom methods new
#' @method group_by GenomicRanges
#' @name group_by-ranges
#' @rdname group_by-ranges
#' @export
#' @examples
#' set.seed(100)
#' df <- data.frame(start = 1:10,
#'                  width = 5,
#'                  gc = runif(10),
#'                  cat = sample(letters[1:2], 10, replace = TRUE))
#' rng <- as_iranges(df)
#' rng_by_cat <- rng %>% group_by(cat)
#' # grouping does not change appearance or shape of Ranges
#' rng_by_cat
#' # a list of symbols
#' groups(rng_by_cat)
#' # ungroup removes any grouping
#' ungroup(rng_by_cat)
#' # group_by works best with other verbs
#' grng <- as_granges(df,
#'                    seqnames = "chr1",
#'                    strand = sample(c("+", "-"), size = 10, replace = TRUE))
#'
#' grng_by_strand <- grng %>% group_by(strand)
#' grng_by_strand
#' # grouping with other verbs
#' grng_by_strand %>% summarise(gc = mean(gc))
#' grng_by_strand %>% filter(gc == min(gc))

#' grng_by_strand %>%
#'   ungroup() %>%
#'   summarise(gc = mean(gc))
#'
#'
group_by.GenomicRanges <- function(.data, ...) {
  new_grouped_gr(.data, ...)
}

#' @method group_by IntegerRanges
#' @export
group_by.IntegerRanges <- function(.data, ...) {
  new_grouped_ir(.data, ...)
}


#' @rdname group_by-ranges
#' @importFrom dplyr ungroup
#' @method ungroup GroupedGenomicRanges
#' @export
ungroup.GroupedGenomicRanges <- function(x, ...) {
  groups <- as.character(unlist(groups(x)))
  capture_groups <- quos(...)
  ungroups <- lapply(capture_groups, function(x) quo_name(x))
  if (length(ungroups) == 0L) {
    return(x@delegate)
  } else {
    groups_update <- setdiff(groups, ungroups)
    groups_update <- syms(groups_update)
    groupings <- make_group_inx(x@delegate, UQS(groups_update))
    new(class(x), delegate = x@delegate, groups = groupings$groups, inx = groupings$inx)
  }
}

#' @method ungroup GroupedIntegerRanges
#' @export
ungroup.GroupedIntegerRanges <- ungroup.GroupedGenomicRanges

#' @method ungroup Ranges
#' @export
ungroup.Ranges <- function(x) x

#' @importFrom dplyr groups
#' @method groups GroupedGenomicRanges
#' @rdname group_by-ranges
#' @export
groups.GroupedGenomicRanges <- function(x) { x@groups }

#' @importFrom dplyr group_vars
#' @method group_vars GroupedGenomicRanges
#' @export
group_vars.GroupedGenomicRanges <- function(x) { as.character(unlist(x@groups)) }

#' @method groups GroupedIntegerRanges
#' @rdname group_by-ranges
#' @export
groups.GroupedIntegerRanges <- groups.GroupedGenomicRanges

#' @method group_vars GroupedIntegerRanges
#' @export
group_vars.GroupedIntegerRanges <- group_vars.GroupedGenomicRanges

#' @method groups Ranges
#' @export
groups.Ranges <- function(x) { NULL }

#' @method group_vars Ranges
#' @export
group_vars.Ranges <- function(x) character(0)