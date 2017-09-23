# ranges-eval-utils.R
# some helpers for 'tidy' NSE on ranges
overscope_ranges <- function(x, envir = parent.frame()) {
  x_env <- as.env(x, envir)
  new_overscope(x_env, top = parent.env(x_env))
}
