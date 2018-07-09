#' @export
tidyr::unnest

#' Expand list-columns in a Ranges object
#'
#' @param data A Ranges object
#' @param ... list-column names to unnest
#' @param .drop Determines whether other list columns will be dropped.
#' By default `unnest` will keep other list columns even if they are nested.
#' @param .id A character vector of length equal to number of list columns.
#' If supplied will create new column(s) with name `.id`
#' identifying the index of the list column (default = NULL).
#' @param .sep Combine name of nested Ranges with name of list column seperated
#' by `.sep`, currently not implemented.
#'
#' @importFrom tidyr unnest
#' @importFrom S4Vectors expand
#' @method unnest GenomicRanges
#'
#' @return a GRanges object with expanded list columns
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


  dots <- rlang::enquos(...)

  dot_names <- unlist(Map(function(.) quo_name(.), dots))

  if (any(!(dot_names %in% colnames(mcols(data))))) {
    stop(paste("Input column(s):",
               paste0(dot_names, collapse = ","), "not found."),
         call. = FALSE)
  }

  if (length(dot_names) == 0L) {
    which_unnest <- names(list_cols)
  } else {
    which_unnest <- intersect(dot_names, names(list_cols))

  }
  if (length(which_unnest) == 0L) {
    stop(paste("Input column(s):",
               paste0(dot_names, collapse = ","), "are not list columns."),
         call. = FALSE)
  }

  if (.drop) {
    list_cols_to_drop <- setdiff(names(list_cols), which_unnest)
    if (length(list_cols_to_drop) > 0) {
      mcols(data) <- mcols(data)[!(names(mcols(data)) %in% list_cols_to_drop)]
    }
  }

  if (!is.null(.id)) {
    if (length(.id) != length(which_unnest)) {
      stop("`.id` does not have same length as number of list columns.",
          call. = FALSE)
    }

    values <- Map(function(col) {
      current <- mcols(data)[[col]]
      v <- names(current)
      if (length(v) == 0) v <- seq_along(current)
      l <- lengths(current)
      inx <- IRanges::IntegerList(Map(function(i) rep(i, l[i]), seq_along(l)))
      IRanges::extractList(v, inx)
      }, which_unnest)


    names(values) <- .id
    values <- do.call("DataFrame", values)
    # expand the ranges first then expand columns
    expand_rng <- S4Vectors::expand(data, which_unnest)
    expand_values <- S4Vectors::expand(values, .id)
    # for some reason if values is single column DF, expand
    # returns a List, with expanded values in last element
    if (length(.id) == 1) {
      expand_values <- expand_values[[.id]]
      expand_values <- DataFrame(expand_values)
      names(expand_values) <- .id
    }
    mcols(expand_rng) <- cbind(mcols(expand_rng), expand_values)
    return(expand_rng)
  }

  return(S4Vectors::expand(data, which_unnest))
}
