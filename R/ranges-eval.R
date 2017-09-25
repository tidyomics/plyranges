# ranges-eval-utils.R
# some helpers for 'tidy' NSE on ranges
overscope_ranges <- function(x, envir = parent.frame()) {
  x_env <- as.env(x, envir)
  new_overscope(x_env, top = parent.env(x_env))
}

ranges_vars <- function(x) {
  x_env <- as.env(x, parent.frame())
  vars_rng <-ls(x_env)
  vars_rng <- vars_rng[!(vars_rng %in% "names")]
  vars_mcols <- ls(parent.env(x_env))
  c(vars_rng, vars_mcols)
}

# dplyr's join syntax uses a function called tbl_vars to get
# variable names, using this function will enable a Ranges to be copied through
# as a data.frame in a join.
tbl_vars.GenomicRanges <- function(x) {
  ranges_vars(x)
}
