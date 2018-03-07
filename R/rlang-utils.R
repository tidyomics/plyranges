# rlang-utils

#' @export
rlang::enquo

#' @export
rlang::quos

#' @export
rlang::UQS

#' @export
rlang::UQ

is_empty_quos <- function(quos) {
  length(quos) == 0L
}
