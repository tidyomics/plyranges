#' Join data by metadata columns
#' 
#' Utilizing DFplyr which enables dplyr verbs for DataFrame
#' objects, these functions enable joins based on metadata
#' of `x` and data in `y`.
#' 
#' @param x Object representing ranges, with metadata
#' columns containing variables for matching
#' @param y A table of data, DataFrame, data.frame or tibble
#' @param ... arguments passed to `left_join`
#' 
#' @return Object representing ranges, with new metadata columns
#'
#' @examples 
#'
#' x <- as_granges(data.frame(
#'   seqnames = "chr1", start = 1:5, width=10,
#'   id=c(1:4,4), foo=letters[1:5]
#' ))
#' y <- DataFrame(id=c(1:2,2,4), bar=letters[6:9])
#'  
#' # metadata join
#' # join_mcols_left(x, y, by="id")
#' 
#' @importFrom dplyr left_join
#' @importFrom rlang .data
#' 
#' @export
join_mcols_left <- function(x, y, ...) {
  x_id <- x
  x_id$.id = factor(seq_along(x))
  new_mcols <- left_join(mcols(x_id), y, ...)
  new_mcols <- arrange(new_mcols, .data[[".id"]])
  new_x <- x[as.integer(new_mcols$.id)]
  new_mcols$.id <- NULL
  mcols(new_x) <- NULL
  mcols(new_x) <- new_mcols
  return(new_x)
}
