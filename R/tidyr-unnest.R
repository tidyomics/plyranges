#' @export
tidyr::unnest

#' Expand list-columns in a Ranges object
#'
#' @param data A Ranges object
#' @param ... list-column names to unnest
#' @param .drop Determines whether other list columns will be dropped.
#' By default \code{unnest} will keep other list columns even if they are nested.
#' @param .id If supplied will create a new column with name \code{.id},
#' identifying the index of the list column.
#' @param .sep Combine name of nested Ranges with name of list column seperated
#' by \code{.sep}/
#'
#' @importFrom tidyr unnest
#' @importFrom S4Vectors expand
#' @method unnest GenomicRanges
#'
#' @examples
#' grng <- as_granges(data.frame(seqnames = "chr1", start = 20:23, width = 1000))
#' grng <- mutate(grng, exon_id = IntegerList(a = 1, b = c(4,5), c = 3, d = c(2,5)))
#' unnest(grng)
#' unnest(grng, .id = "name")
#' @rdname ranges-unnest
#' @export
unnest.GenomicRanges <- function(data, ..., .drop = FALSE, .id = NULL, .sep = NULL) {

  list_cols_pos <- unlist(Map(function(.) is(., "List"), mcols(data)))

  list_cols <- Filter(isTRUE, list_cols_pos)

  if (length(list_cols) == 0L) {
    stop("No list columns to unnest.", call. = FALSE)
  }


  dots <- quos(...)

  dot_names <- unlist(Map(function(.) quo_name(.), dots))


  if (length(dot_names) == 0L) {
    which_unnest <- names(list_cols)
  } else {
    which_unnest <- intersect(dot_names, names(list_cols))

  }
  if (length(which_unnest) == 0L) {
    stop(paste("Input column(s):",
               paste0(dot_names, collapse = ","), "not found"), call. = FALSE)
  }
  if (.drop) {
    list_cols_to_drop <- setdiff(names(list_cols), which_unnest)
    if (length(list_cols_to_drop) > 0) {
      mcols(data) <- mcols(data)[!(names(mcols(data)) %in% list_cols_to_drop)]
    }
  }

  unnest_rng <- S4Vectors::expand(data, which_unnest)
  if (!is.null(.id)) {
    if (length(which_unnest) == 1L) {
      values <- names(mcols(data)[[which_unnest]])
      if (is.null(values)) {
        values <- seq_along(mcols(data)[[which_unnest]])
      }
      lengths <- lengths(mcols(data)[[which_unnest]])
      mcols(unnest_rng)[[.id]] <- Rle(values, lengths)
    }
  }
  unnest_rng
}
