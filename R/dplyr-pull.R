#' Extract a single column from a Ranges object as a vector
#'
#' @param .data a `Ranges` object
#' @param var A variable specified as:
#'   * a literal variable name
#'   * a positive integer, giving the position counting from the left. In this
#'   case order is start, end, width, (strand, seqnames), gc and score. 
#'   * a negative integer, giving the position counting from the right.
#'   The default returns the last column (on the assumption that's the column 
#'   you've created most recently). This argument is taken by expression and 
#'   supports quasiquotation (you can unquote column names and column locations).
#' @param name An optional parameter that specifies the column to be used as 
#'   names for a named vector. Specified in a similar manner as `var`.
#' @param ... For use by methods.
#'
#'
#' @importFrom dplyr pull
#' @importFrom rlang enquo quo_is_missing
#' @importFrom tidyselect eval_select
#' @seealso [dplyr::pull()]
#' @name pull-ranges
#' @rdname pull-ranges
#'
#' @examples
#' df <- data.frame(start = 1:10,
#'                  width = 5,
#'                  seqnames = "seq1",
#'                  strand = sample(c("+", "-", "*"), 10, replace = TRUE),
#'                  gc = runif(10),
#'                  score = rpois(10, 2))
#' rng <- as_granges(df)
#'
#' # Pull parts of the range
#' pull(rng, start)
#' # equivalent to start(rng)
#'
#' # Pull by column name
#' pull(rng, gc)
#' pull(rng, score)
#'
#' # Pull by position (positive from left, negative from right)
#' pull(rng, 1)    # First metadata column
#' pull(rng, -1)   # Last metadata column (default)
#' pull(rng, -2)   # Second to last metadata column
#'
#' # Pull with names
#' pull(rng, score, name = gc)
#'
#' @method pull Ranges
#' @exportS3Method dplyr::pull
pull.Ranges <- function(.data, var = -1, name = NULL, ...) {
  
  var <- tidyselect::vars_pull(tbl_vars(.data), !!rlang::enquo(var))
  name <- rlang::enquo(name)
  env <- overscope_ranges(.data) 
  
  pulled <- eval_tidy(rlang::ensym(var), env)
  if (rlang::quo_is_null(name)) {
    return(pulled)
  }
  
  name <- tidyselect::vars_pull(tbl_vars(.data), !!name)
  named <- eval_tidy(rlang::ensym(name), env)
  
  rlang::set_names(pulled, nm = named)
  
}
  

#' @method pull DelegatingGenomicRanges
#' @export
pull.DelegatingGenomicRanges <- function(.data, var = -1, name = NULL, ...) {
  pull(.data@delegate, var = {{var}}, name = {{name}}, ...)
}

#' @method pull DelegatingIntegerRanges
#' @export
pull.DelegatingIntegerRanges <- pull.DelegatingGenomicRanges

#' @method pull GroupedGenomicRanges
#' @export
pull.GroupedGenomicRanges <- function(.data, var = -1, name = NULL, ...) {
  pull(.data@delegate, var = {{var}}, name = {{name}}, ...)
}

#' @method pull GroupedIntegerRanges
#' @export
pull.GroupedIntegerRanges <- pull.GroupedGenomicRanges