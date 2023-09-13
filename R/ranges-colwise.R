# Various column wise helpers for ranges

#' Tools for working with named Ranges
#'
#' @param .data a Ranges object
#' @param var Name of column to use for names
#' 
#' @details The function `names_to_column()` and `id_to_column()` always places
#' `var` as the first column in `mcols(.data)`, shifting all other columns 
#' to the left. The `id_to_column()` creates a column with sequential row
#' identifiers starting at 1, it will also remove any existing names. 
#' 
#' @return 
#' Returns a Ranges object with empty names
#' 
#' @examples 
#' ir <- IRanges::IRanges(start = 1:3, width = 4, names = c("a", "b", "c"))
#' remove_names(ir)
#' ir_noname <- names_to_column(ir)
#' ir_noname
#' ir_with_id <- id_to_column(ir)
#' ir_with_id
#' @rdname ranges-names
#' @export
remove_names <- function(.data) {
  stopifnot(is(.data, "Ranges"))
  if (is.null(names(.data))) return(.data)
  names(.data) <- NULL
  .data
}


#' @rdname ranges-names
#' @export
names_to_column <- function(.data, var = "name") {
  stopifnot(is(.data, "Ranges"))
  if (var %in% names(mcols(.data))) {
    stop(paste("Column", var, "already exists in .data."))
  }
  if (is.null(names(.data))) {
    stop("No names present in .data")
  }
  r_names <- names(.data)
  .data <- mutate(.data, `:=`(!!rlang::sym(var), r_names))
  .data <- select(.data, !!var, tidyselect::everything())
  remove_names(.data)
}

#' @rdname ranges-names
#' @export
id_to_column <- function(.data, var = "id") {
  stopifnot(is(.data, "Ranges"))
  if (var %in% names(mcols(.data))) {
    stop(paste("Column", var, "already exists in .data."))
  }
  r_id <- seq_len(length(.data))
  .data <- mutate(.data, `:=`(!!var, r_id))
  .data <- select(.data, !!var, tidyselect::everything())
  remove_names(.data)
}