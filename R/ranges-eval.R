# ranges-eval-utils.R
# some helpers for 'tidy' NSE on ranges

are_bioc_pkgs_scoped <- function(pkgs = c("BiocGenerics", "IRanges", "S4Vectors")) {
  all(vapply(paste0("package:", pkgs), 
             rlang::is_scoped, 
             logical(1)))
}

#' @importFrom methods getGeneric getGenerics
bioc_generics <- function(pkgs = c("BiocGenerics", "IRanges", "S4Vectors")) {
  pkgs <- lapply(pkgs, asNamespace)
  generics <-  getGenerics(pkgs)
  fn <- mapply(getGeneric, 
               f = generics@.Data, 
               package = generics@package, 
               SIMPLIFY = FALSE)
  Filter(function(x) {
    fun = try(x@default, silent = TRUE)
    if (is(fun, "try-error")) FALSE
    !is.primitive(fun)
  },
  fn)
}

scope_plyranges <- function(env, generics) {
  tail <- env
  nms <- character(0)
  # recurse through parent environments until we get to the empty env
  while(!identical(tail, rlang::empty_env())) {
    env_nms <- rlang::env_names(tail)
    nms <- unique(c(nms, intersect(names(generics), env_nms)))
    tail <- rlang::env_parent(tail)
  }
  nms <- setdiff(nms, rlang::env_names(rlang::global_env()))
  generics <- generics[nms]
  
  child <- rlang::child_env(env,
                            UQS(generics))
  child
}

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
  generics <- bioc_generics()
  for (i in seq_along(update)) {
    quo <- dots[[i]]
    qenv <- rlang::quo_get_env(quo)
    quo <- rlang::quo_set_env(quo, scope_plyranges(qenv, generics))
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
  if (rlang::env_has(parent_env, "start")) {
    .data <- rlang::env_get(parent_env, "start")
    if (is(.data, "IntegerList")) {
      return(BiocGenerics::lengths(.data))
    } else {
      return(length(.data))
    }
  }
  stop("This function should not be called directly")
}

# dplyr's join syntax uses a function called tbl_vars to get
# variable names, using this function will enable a Ranges to be copied through
# as a data.frame in a join.
tbl_vars.GenomicRanges <- function(x) {
  ranges_vars(x)
}


detach_depends <- function() {
  sapply(paste0("package:", 
                c("GenomicRanges",
                  "GenomeInfoDb",
                  "IRanges",
                  "S4Vectors",
                  "BiocGenerics")),
         detach, 
         character.only = TRUE)
}

reattach_depends <- function() {
  sapply(c("GenomicRanges",
           "GenomeInfoDb",
           "IRanges",
           "S4Vectors",
           "BiocGenerics"),
         library,
         character.only = TRUE,
         quietly = TRUE)
}
