#' Group a Ranges by one or more variables
#'
#' @param .data a Ranges object.
#' @param ... Variable names to group by. These can be either metadata columns
#' or the core variables of a Ranges.
#' @param add if `.data` is already a GroupedRanges object, when add = FALSE 
#' the (default), `group_by()` will override existing groups. If add = TRUE, 
#' additional groups will be added.
#' @param x a GroupedRanges object.
#'
#' @description The function `group_by` takes a Ranges object and defines
#' groups by one or more variables. Operations are then performed on the Ranges
#' by their "group". `ungroup()` removes grouping.
#'
#' @details
#' `group_by()` creates a new object of class `GroupedGenomicRanges` if
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
group_by.GenomicRanges <- function(.data, ..., add = FALSE) {
  new_grouping(.data, ...)
}

#' @method group_by IntegerRanges
#' @export
group_by.IntegerRanges <- function(.data, ..., add = FALSE) {
  new_grouping(.data, ..., target = "GroupedIntegerRanges")
}

#' @method group_by GroupedGenomicRanges
#' @export
group_by.GroupedGenomicRanges <- function(.data, ..., add = FALSE) {
  if (add) {
    
  }
  
  new_grouping(.data@delegate, ..., target = "GroupedIntegerRanges")
}

#' @rdname group_by-ranges
#' @importFrom dplyr ungroup
#' @method ungroup GroupedGenomicRanges
#' @export
ungroup.GroupedGenomicRanges <- function(x, ...) {
  ungroups <- enquos(...)
  ungroups <- rlang::quos_auto_name(ungroups)
  if (length(ungroups) == 0L) {
    return(x@delegate)
  } else {
    ungroups <- lapply(ungroups, function(.) rlang::quo(-!!.))
    groups_update <- tidyselect::vars_select(group_vars(x), !!!ungroups)
    if (length(groups_update) == 0) {
      return(x@delegate)
    }
    
    groups_update <- syms(groups_update)
    new_grouping(x@delegate, !!!groups_update)
  }
}

#' @method ungroup GroupedIntegerRanges
#' @export
ungroup.GroupedIntegerRanges <- ungroup.GroupedGenomicRanges

#' @method ungroup Ranges
#' @export
ungroup.Ranges <- function(x, ...) x

#' @importFrom dplyr groups
#' @method groups GroupedGenomicRanges
#' @rdname group_by-ranges
#' @export
groups.GroupedGenomicRanges <- function(x)  syms(colnames(x@group_keys))

#' @importFrom dplyr group_vars
#' @method group_vars GroupedGenomicRanges
#' @export
group_vars.GroupedGenomicRanges <- function(x) colnames(x@group_keys) 

#' @method groups GroupedIntegerRanges
#' @rdname group_by-ranges
#' @export
groups.GroupedIntegerRanges <- groups.GroupedGenomicRanges

#' @method group_vars GroupedIntegerRanges
#' @export
group_vars.GroupedIntegerRanges <- group_vars.GroupedGenomicRanges

#' @method groups Ranges
#' @export
groups.Ranges <- function(x)  NULL 

#' @method group_vars Ranges
#' @export
group_vars.Ranges <- function(x) character(0)

#' @method group_keys GroupedGenomicRanges
#' @export
#' @importFrom dplyr group_keys
group_keys.GroupedGenomicRanges <- function(.data, ...) {
  .data@group_keys
}

#' @method group_keys GroupedIntegerRanges
#' @export
group_keys.GroupedIntegerRanges <- group_keys.GroupedGenomicRanges

#' @method group_keys Ranges
#' @export
group_keys.Ranges <- function(.data, ...) {
  if (length(enquos(...)) == 0) {
    return(new("DFrame", nrows = 1L))
  }
  NextMethod(group_by(.data, ...))
}

#' @method group_indices GroupedGenomicRanges
#' @export
#' @importFrom dplyr group_indices
group_indices.GroupedGenomicRanges <- function(.data, ...) {
  .data@group_indices
} 

#' @method group_indices GroupedIntegerRanges
#' @export
group_indices.GroupedIntegerRanges <- group_indices.GroupedGenomicRanges


#' @method group_indices Ranges
#' @export
group_indices.Ranges <- function(.data, ...) {
  if (length(enquos(...)) == 0) {
    return(rep.int(1, length(.data)))
  }
  NextMethod(group_by(.data, ...))
}


.group_rows <- function(.data) {
  as(unname(S4Vectors::split(
    seq_len(length(.data@delegate)),
    .data@group_indices
  )),
  "IntegerList"
  )
}

#' @method group_data GroupedGenomicRanges
#' @export
#' @importFrom dplyr group_data  
group_data.GroupedGenomicRanges <- function(.data) {
  S4Vectors::DataFrame(
    .data@group_keys,
    .rows = .group_rows(.data))
}

#' @method group_data GroupedIntegerRanges
#' @export
group_data.GroupedIntegerRanges <- group_data.GroupedGenomicRanges

#' @method group_data Ranges
group_data.Ranges <- function(.data) {
  .rows <- as(seq_len(length(.data)), "IntegerList")
  S4Vectors::DataFrame(.rows = .rows)
}

#' @method group_split GroupedGenomicRanges
#' @export
#' @importFrom dplyr group_split
group_split.GroupedGenomicRanges <- function(.data, ..., keep = TRUE) {
  if (length(enquos(...)) > 0) {
    warning("Ignoring arguments to `...` 
            and using existing group structure")
  }
  
  rng <- .data@delegate 
  
  if (!keep) {
    vars_drop <- lapply(group_vars(.data), function(.) rlang::quo(-!!.))
    rng <- select(rng, !!!vars_drop)
  } 
  
  unname(S4Vectors::split(rng, .data@group_indices))
}

#' @method group_split GroupedIntegerRanges
#' @export
group_split.GroupedIntegerRanges <- group_split.GroupedGenomicRanges

#' @method group_split Ranges
#' @export
group_split.Ranges <- function(.data, ..., keep = TRUE) {
  if (length(enquos(...)) == 0) {
    return(as(.data, "List"))
  }
  NextMethod(group_by(.data, ...))
}