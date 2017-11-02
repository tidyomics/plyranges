# rlang-utils

#' @export
rlang::enquo

#' @export
rlang::quos

#' @export
rlang::UQS

#' @export
rlang::UQ

#' @export
rlang::overscope_clean

#' @export
rlang::overscope_eval_next

#' @export
rlang::new_overscope

#' @export
rlang::eval_bare

#' @export
rlang::eval_tidy

#' @export
rlang::syms

#' @export
rlang::`:=`

#' @export
rlang::env_bind

is_empty_quos <- function(quos) {
  length(quos) == 0L
}

#' @export
tibble::as_tibble
