# ranges-eval-utils.R
# some helpers for 'tidy' evaluation on ranges


#' Create an overscoped environment from a Ranges object
#' 
#' @param x a Ranges object
#' @param envir the environment to place the Ranges in (default = `parent.frame()`)
#' 
#' @details This is the backend for non-standard evaluation in `plyranges`.
#' 
#' @seealso [rlang::new_data_mask()], [rlang::eval_tidy()]
#' @return an environment
#' 
#' @export
overscope_ranges <- function(x, envir = parent.frame()) {
  UseMethod("overscope_ranges")
}

overscope_ranges.Ranges <- function(x, envir = parent.frame()) {
  env <- as.env(x, envir)
  
  new_data_mask(env, top = parent.env(env))
}

overscope_ranges.DelegatingGenomicRanges <- function(x, envir = parent.frame()) {
  overscope_ranges(x@delegate, envir)
}

overscope_ranges.DelegatingIntegerRanges <- overscope_ranges.DelegatingGenomicRanges

overscope_ranges.GroupedGenomicRanges <- function(x, envir = parent.frame()) {
  env <- as.env(x@delegate, 
                  envir, 
                  tform = function(col) unname(IRanges::extractList(col, x@inx)))
  new_data_mask(env, top = parent.env(env))
}

overscope_ranges.GroupedIntegerRanges <- overscope_ranges.GroupedGenomicRanges



#' @importFrom rlang env_bind := new_data_mask eval_tidy
overscope_eval_update <- function(overscope, dots, bind_envir = TRUE) {
  update <- vector("list", length(dots))
  names(update) <- names(dots)
  for (i in seq_along(update)) {
    quo <- dots[[i]]
    update[[i]] <- eval_tidy(quo, data = overscope)
    # sometimes we want to compute on previously constructed columns
    # we can do this by binding the evaluated expression to
    # the overscope environment
    if (bind_envir) {
      new_col <- names(dots)[[i]]
      rlang::env_bind(overscope, !!new_col := update[[i]])
    }
  }
  return(update)
}


# Port of dplyrs `n` function
# It works by searching for a vector in the overscope environment
# and calling length on it.

#' Compute the number of ranges in each group.
#'
#' @description This function should only be used
#' within `summarise()`, `mutate()` and `filter()`.
#'
#' @examples
#' ir <- as_iranges(
#'                  data.frame(start = 1:10,
#'                             width = 5,
#'                             name = c(rep("a", 5), rep("b", 3), rep("c", 2))
#'                             )
#'                 )
#' by_names <- group_by(ir, name)
#' summarise(by_names, n = n())
#' mutate(by_names, n = n())
#' filter(by_names, n() >= 3)
#' @return `n()` will only be evaluated inside a function call, where it
#' returns an integer.
#'
#' @importFrom rlang env_get env_parent
#' @export
n <- function() {
  up_env <- parent.frame()
  parent_env <- rlang::env_parent(up_env)
  if (rlang::env_has(parent_env, "start")) {
    .data <- rlang::env_get(parent_env, "start")
    if (is(.data, "IntegerList")) {
      return(lengths(.data))
    } else {
      return(length(.data))
    }
  }
  stop("This function should not be called directly")
}


#' Compute the number of distinct unique values in a vector or List
#' 
#' @param var a vector of values
#' @return an integer vector 
#' 
#' @description This is a wrapper to `length(unique(x))` or 
#' `lengths(unique(x))` if `x` is a List object
#' 
#' @examples 
#' x <- CharacterList(c("a", "b", "c", "a"),  "d")
#' n_distinct(x)
#' n_distinct(unlist(x))
#' 
#' @export
n_distinct <- function(var) {
  if (inherits(var, "List")) {
    return(lengths(unique(var)))
  } else {
    return(length(unique(var)))
  }
}


is_empty_quos <- function(quos) {
  length(quos) == 0L
}


#' @importFrom dplyr tbl_vars
tbl_vars.Ranges <- function(x) {
  c("start", "end", "width", names(mcols(x)))
}

tbl_vars.GenomicRanges <- function(x) {
  c("start", "end", "width", "strand", "seqnames", names(mcols(x)))
}