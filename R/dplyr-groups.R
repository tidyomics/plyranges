#' Group a Ranges by one or more variables
#'
#' @param .data a Ranges object.
#' @param ... Variable names to group by. These can be either metadata columns
#' or the core variables of a Ranges.
#' @param x a GroupedRanges object.
#'
#' @description The function \code{group_by} takes a Ranges object and defines
#' groups by one or more variables. Operations are then performed on the Ranges
#' by their "group". \code{ungroup()} removes grouping.
#'
#' @details
#' \code{group_by()} creates a new object of class \code{GRangesGrouped} if
#' the input is a \code{GRanges} object or an object of class \code{IRangesGrouped}
#' if the input is a \code{IRanges} object. Both of these classes contain a slot
#' called \code{groups} corresponding to the names of grouping variables. They
#' also inherit from their parent classes, \code{Ranges} and \code{GenomicRanges}
#' respectively. \code{ungroup()} removes the grouping and will return
#' either a \code{GRanges} or \code{IRanges} object.
#'
#' @section Accessors:
#' To return grouping variables on a grouped Ranges use either
#' \itemize{
#'   \item{\code{groups(x)}}{Returns a list of symbols}
#'   \item{\code{group_vars(x)}}{Returns a character vector}
#' }
#'
#' @return The \code{group_by()} function will return a GroupedRanges object.
#' These have the same appearance as a regular Ranges object but with an
#' additional groups slot.
#'
#'
#' @importFrom dplyr group_by
#' @importFrom rlang quo_name quos syms
#' @importFrom methods new
#' @method group_by GRanges
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
group_by.GRanges <- function(.data, ...) {
  capture_groups <- quos(...)
  groups <- lapply(capture_groups, function(x) quo_name(x))
  groups <- syms(groups)
  new("GRangesGrouped", .data, groups = groups)

}



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
    rng <- IRanges(x)
    if (!is.null(mcols(x))) {
      mcols(rng) <- mcols(x)
    }
    return(rng)
  } else {
    groups_update <- setdiff(groups, ungroups)
    groups_update <- syms(groups_update)
    new("IRangesGrouped", x, groups = groups_update)
  }
}


#' @importFrom dplyr groups
#' @method groups GRangesGrouped
#' @rdname group_by-ranges
#' @export
groups.GRangesGrouped <- function(x) { x@groups }

#' @importFrom dplyr group_vars
#' @method group_vars GRangesGrouped
#' @rdname group_by-ranges
#' @export
group_vars.GRangesGrouped <- function(x) { as.character(unlist(x@groups)) }

#' @method groups IRangesGrouped
#' @rdname group_by-ranges
#' @export
groups.IRangesGrouped <- function(x) { x@groups }

#' @method group_vars IRangesGrouped
#' @rdname group_by-ranges
#' @export
group_vars.IRangesGrouped <- function(x) { as.character(unlist(x@groups)) }
