

#' Expand list-columns in a Ranges object
#'
#' @param data A Ranges object
#' 
#' @param ... list-column names to expand then unlist
#' @param .drop Should additional list columns be dropped (default = FALSE)?
#' By default `expand_ranges()` will keep other list columns even if they are nested.
#' @param .id A character vector of length equal to number of list columns.
#' If supplied will create new column(s) with name `.id`
#' identifying the index of the list column (default = NULL).
#' @param .keep_empty If a list-like column contains empty elements, should
#' those elements be kept? (default = FALSE)
#' @param .recursive If there are multiple list-columns, should the columns be 
#' treated as parallel? If FALSE each column will be unnested recursively,
#' otherwise they are treated as parallel, that is each list column has
#' identical lengths. (deafualt = FALSE)
#'
#' @importFrom S4Vectors expand
#'
#' @return a GRanges object with expanded list columns
#' @examples
#' grng <- as_granges(data.frame(seqnames = "chr1", start = 20:23, width = 1000))
#' grng <- mutate(grng, 
#'                exon_id = IntegerList(a = 1, b = c(4,5), c = 3, d = c(2,5))
#'                )
#' expand_ranges(grng)
#' expand_ranges(grng, .id = "name")
#' 
#' # empty list elements are not preserved by default
#' grng <- mutate(grng, 
#'                exon_id = IntegerList(a = NULL, b = c(4,5), c= 3, d = c(2,5))
#'                )
#' expand_ranges(grng)
#' expand_ranges(grng, .keep_empty = TRUE)
#' expand_ranges(grng, .id = "name", .keep_empty = TRUE)
#' 
#' @rdname ranges-expand
#' @export
expand_ranges <- function(data, ..., .drop = FALSE, .id = NULL, .keep_empty = FALSE, .recursive = FALSE) {
  
  list_cols <- get_list_cols(data)
  which_unnest <- unnest_cols(data, list_cols, ...)
  
  if (.drop) {
    list_cols_to_drop <- setdiff(names(list_cols), which_unnest)
    if (length(list_cols_to_drop) > 0) {
      mcols(data) <- mcols(data)[!(names(mcols(data)) %in% list_cols_to_drop)]
    }
  }
  
  if (ncol(mcols(data)) == length(which_unnest)) {
    expand_rng <- S4Vectors::expand(data, 
                                    keepEmptyRows = .keep_empty,
                                    recursive = .recursive)
  } else {
    expand_rng <- S4Vectors::expand(data, 
                                    which_unnest, 
                                    keepEmptyRows = .keep_empty,
                                    recursive = .recursive)
  }
  
  
  
  if (!is.null(.id)) {
    if (length(.id) != length(which_unnest)) {
      stop("`.id` does not have same length as number of list columns.",
           call. = FALSE)
    }
    
    values <- lapply(mcols(data)[, which_unnest, drop = FALSE],
                     extract_index,
                     .keep_empty = .keep_empty)
    
    names(values) <- .id
    values <- do.call("DataFrame", values)
    # for some reason if values is single column DF, expand
    # returns a List, with expanded values in last element
    if (length(.id) == 1) {
      expand_values <- expand(values, 
                              .id, 
                              keepEmptyRows = .keep_empty,
                              recursive = .recursive)
      expand_values <- expand_values[[.id]]
      expand_values <- DataFrame(expand_values)
      names(expand_values) <- .id
    } else {
      expand_values <- S4Vectors::expand(values, 
                                         keepEmptyRows = .keep_empty,
                                         recursive = .recursive )
    }
    
    mcols(expand_rng) <- cbind(mcols(expand_rng), expand_values)
  }
  
  return(expand_rng)
  
}


get_list_cols <- function(data) {
  list_cols_pos <- unlist(Map(function(.) is(., "List"), mcols(data)))
  
  list_cols <- Filter(isTRUE, list_cols_pos)
  
  if (length(list_cols) == 0L) {
    stop("No list columns to unnest.", call. = FALSE)
  }
  
  list_cols
}

unnest_cols <- function(data, list_cols, ...) {
  
  dots <- set_dots_unnamed(...)
  
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
  
  return(which_unnest)
}

extract_index <- function(col, .keep_empty) {

  v <- names(col)
  if (length(v) == 0) v <- seq_along(col)
  l <- lengths(col)
  if (.keep_empty) l[l == 0] <- 1
  inx <- IRanges::IntegerList(Map(function(i) rep(i, l[i]), seq_along(l)))
  IRanges::extractList(v, inx)
  
}
