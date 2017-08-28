#' Subset a GRanges by row
#'
#' @param .data a \code{GRanges} object
#' @param row_values integer row values
#' @param coverage slice based on coverage values? default = FALSE
#'
#' @return a \code{GRanges} object
#' this conflicts with IRanges::slice, and perhaps should do something else
slice.GRanges <- function(.data, row_values, coverage = FALSE) {

  row_values <- enquo(row_values)
  row_index <- quo_name(row_values)
  # little hacky to get n() equivalent
  row_index <- gsub("n\\(\\)", length(.data), row_index)
  is_slice <- grepl(":", row_index)
  if (is_slice) {
    start_row <- as.integer(sub("^(.*?):(.*)", "\\1", row_index))
    end_row <- sub("^(.*?):(.*)", "\\2", row_index)
    final_range <- seq.int(as.integer(start_row), as.integer(end_row))
  } else {
    final_range <- as.integer(row_index)
  }

  if (coverage) {
    message("not quite sure yet")
  }

  .data[final_range, ]
}
