#' Select metadata variables by name
#'
#' @param .data a \code{GRanges} object
#' @param ... One or more column names. Note that these act on the
#' metadata slot of the GRanges object and preserve the
#' structure of the GRanges object, hence the names strand, start, end,
#' width are invalid.
#' # @param drop drop the GRanges column and return the metadata columns only?
#' @return a GRanges object
#' @seealso \link[dplyr]{select}
select.GRanges <- function(.data, ...) {
  dots <- quos(...)

  col_names <- unlist(lapply(dots, quo_name))

  invalid_cols <- any(grepl("strand|start|end|width", col_names))
  if (invalid_cols) {
    stop("Invalid selection, only metadata columns may be selected.")
  }

  # coerce to a data.frame
  mcols(.data) <- select(as.data.frame(mcols(.data)),
                         UQS(dots))

  .data
  # col_slice <- grep(":", col_names)
  #
  # if (length(col_slice) > 0) {
  #   start_col <- sub("^(.*?):(.*)",
  #                    "\\1",
  #                    col_names[col_slice])
  #   end_col <- sub("^(.*?):(.*)", "\\2", col_names[col_slice])
  #   col_slice_index <- lapply(seq_along(col_slice),
  #                             function(i) {
  #                               start <- grep(start_col[i], mcol_names)
  #                               end <- grep(end_col[i], mcol_names)
  #                               seq.int(start, end, by = 1)
  #                               })
  #   col_slice_index <- unlist(col_slice_index)
  #   col_names <- col_names[-col_slice]
  # } else {
  #   col_slice_index <- integer(0)
  # }
  #
  # if (length(col_names) > 0) {
  #   print(col_names)
  #   col_index <- vapply(col_names,
  #                       function(x)  {
  #                         grep(sub("^-", "", x), mcol_names) *
  #                           ifelse(grepl("^-", x), -1L, 1L)
  #                       },
  #                       integer(1))
  #   col_drops <- sign(col_index) == -1
  # } else {
  #   col_index <- integer(0)
  # }
  #
  # if (any(col_drops)) {
  #   .data <- .data[, col_index[col_drops]]
  #   col_index <- col_index[!col_drops]
  # }
  #
  # .data[, c(col_slice_index,col_index)]

}
