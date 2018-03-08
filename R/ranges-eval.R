# ranges-eval-utils.R
# some helpers for 'tidy' NSE on ranges
overscope_ranges <- function(x, envir = parent.frame()) {
  x_env <- as.env(x, envir)
  os <- new_overscope(x_env, top = parent.env(x_env))
  return(os)
}

#' @importFrom rlang env_bind := new_overscope overscope_eval_next
overscope_eval_update <- function(overscope, dots, bind_envir = TRUE) {
  update <- vector("list", length(dots))

  names(update) <- names(dots)
  for (i in seq_along(update)) {
    update[[i]] <- overscope_eval_next(overscope, dots[[i]])
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

ranges_vars <- function(x) {
  x_env <- as.env(x, parent.frame())
  vars_rng <-ls(x_env)
  vars_rng <- vars_rng[!(vars_rng %in% "names")]
  vars_mcols <- ls(parent.env(x_env))
  c(vars_rng, vars_mcols)
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
  .data <- try(rlang::env_get(parent_env, "start"),
               silent = TRUE)
  if (is(.data, "try-error")) {
    stop("This function should not be called directly.")
  }
  return(length(.data))
}

# dplyr's join syntax uses a function called tbl_vars to get
# variable names, using this function will enable a Ranges to be copied through
# as a data.frame in a join.
tbl_vars.GenomicRanges <- function(x) {
  ranges_vars(x)
}
