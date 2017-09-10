# arithemtic-ranges.R
#' @importFrom rlang quo
stretch <- function(x, extend = 0L) {

  expr <- enquo(x)
  print(expr)
  print(quo_name(expr))
  if (grepl("start", quo_name(expr))) {
    quos(width = width + UQ(extend))
  } else if (grepl("end", quo_name(expr))) {
    quos(start = end - width + UQ(extend) + 1L)
  } else if (grepl("center", quo_name(expr))) {
    midpt <- eval_tidy(expr)
    if (sign(extend) == -1) {
      quos(start = start + 2*midpt, end = start)
    } else {
      quos(start = end, end = end + 2*midpt)
    }
  } else {
    !!!quos(shift(x, !!extend))
  }

}

anchor_start <- function() {
  stop("anchor_start should only be called inside another function",
       .call = FALSE)
}

anchor_end <- function(x) {
  end(x)
}

anchor_center <- function(x) {
  (end(x) - start(x)) / 2
}

