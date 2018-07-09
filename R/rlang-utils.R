# rlang-utils
# allows programming with plyranges

#' @export
rlang::enquo

#' @export
rlang::UQS

#' @export
rlang::UQ

is_empty_quos <- function(quos) {
  length(quos) == 0L
}
