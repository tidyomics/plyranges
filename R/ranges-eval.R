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

#' @export
overscope_ranges.Ranges <- function(x, envir = parent.frame()) {
  env <- overscope_env(x, envir)

  new_data_mask(env, top = parent.env(env))
}

# The (non-metadata) parallel-slot names of a Ranges object, i.e.
# seqnames/start/end/width/strand. Introduced because of
# (https://github.com/Bioconductor/S4Vectors/issues/140). Several
# call sites treat parallelVectorNames() as "the core slots" 
# if upstream issue is fixed.
core_vector_names <- function(x) {
  setdiff(S4Vectors::parallelVectorNames(x), names(mcols(x, use.names = FALSE)))
}

#' Build the two-tier data-mask environment for a Ranges object
#'
#' @description
#' Internal backend for [overscope_ranges()]. Constructs the environment used
#' for non-standard evaluation: a child tier binding the fixed parallel slots
#' (seqnames/start/end/width/strand) enclosed by a parent tier binding the
#' metadata columns, mirroring the layout produced by [IRanges::as.env()].
#'
#' @param x a Ranges object.
#' @param envir the enclosing environment for the returned environment.
#' @param tform a function applied to each bound column, used to split columns
#'   into a List for grouped Ranges; defaults to [identity()].
#'
#' @return An environment whose child tier binds the fixed parallel slots and
#'   whose parent (enclosing) tier binds the metadata columns.
#'
#' @seealso [overscope_ranges()], [IRanges::as.env()]
#' @noRd
overscope_env <- function(x, envir, tform = identity) {
  mcols_env <- as.env(mcols(x, use.names = FALSE), envir, tform)
  x_bare <- x
  mcols(x_bare) <- NULL
  env <- as.env(x_bare, envir, tform)
  parent.env(env) <- mcols_env
  env$.. <- x
  env
}

#' @export
overscope_ranges.DelegatingGenomicRanges <- function(x, envir = parent.frame()) {
  overscope_ranges(x@delegate, envir)
}

#' @export
overscope_ranges.DelegatingIntegerRanges <- overscope_ranges.DelegatingGenomicRanges

#' @export
overscope_ranges.GroupedGenomicRanges <- function(x, envir = parent.frame()) {
  env <- overscope_env(x@delegate,
                       envir,
                       tform = function(col) unname(S4Vectors::splitAsList(col, x@group_indices)))
  new_data_mask(env, top = parent.env(env))
}


#' @export
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
#' @importFrom BiocGenerics unique
#' @export
n_distinct <- function(var) {
  if (inherits(var, "List")) {
    return(lengths(BiocGenerics::unique(var)))
  } else {
    return(length(BiocGenerics::unique(var)))
  }
}

is_empty_quos <- function(quos) {
  length(quos) == 0L
}
