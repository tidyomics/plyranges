# ranges-eval-utils.R
# some helpers for 'tidy' NSE on ranges
overscope_ranges <- function(x, envir = parent.frame(), bind_data = FALSE) {
  x_env <- as.env(x, envir)
  os <- new_overscope(x_env, top = parent.env(x_env))
  if (bind_data) {
    rlang::env_bind(os, .data := x)
  }
  return(os)
}

#' @importFrom rlang env_bind :=
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


is_n <- function(dots) {
  qnames <- unlist(lapply(dots, rlang::quo_name))
  any(grepl("n\\(\\)", qnames))
}

#' @importFrom rlang quo_name get_env parse_quosure
check_n <- function(dots) {
  if (rlang::is_quosure(dots)) {
    qnames <- rlang::quo_name(dots)
    modify_n <- gsub("n\\(\\)", "length\\(.data\\)", qnames)
    return(rlang::parse_quosure(modify_n, rlang::get_env(dots)))

  } else {
    qnames <- unlist(lapply(dots, rlang::quo_name))
    which_n <- grep("n\\(\\)", qnames)
    if (length(which_n) == 0L) {
      return(dots)
    } else {
      modify_n <- gsub("n\\(\\)", "length\\(.data\\)", qnames)

      for (i in seq_along(which_n)) {
        update <- which_n[i]
        dots[[update]] <- rlang::parse_quosure(modify_n[update],
                                               env = rlang::get_env(dots[[update]]))
      }
      return(dots)
    }
  }
}

# dplyr's join syntax uses a function called tbl_vars to get
# variable names, using this function will enable a Ranges to be copied through
# as a data.frame in a join.
tbl_vars.GenomicRanges <- function(x) {
  ranges_vars(x)
}
