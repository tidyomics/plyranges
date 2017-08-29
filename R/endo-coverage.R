#' GRanges-friendly coverage method
#'
#' @param x a \code{GRanges} object
#' @param shift shift
#' @param width width
#' @param weight weight
#' @param method method
#'
#' @return an expanded GRanges object with a score column corresponding to
#' the coverage value over that interval. Note that compute_coverage
#' drops metadata associated with the orginal ranges.
#' @seealso \link[IRanges]{coverage}
#' @importFrom IRanges coverage
#' @export
setGeneric("compute_coverage",
           function(x, shift = 0L, width = NULL, weight = 1L,
                    method = c("auto", "sort", "hash")) {
             standardGeneric("compute_coverage")
             })

setMethod("compute_coverage", "GRanges",
          function(x, shift = 0L, width = NULL, weight = 1L,
                   method = c("auto", "sort", "hash")) {
            GRanges(coverage(x, shift, width, weight, method))
          })

