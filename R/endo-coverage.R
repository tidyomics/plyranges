#' Ranges friendly coverage method
#'
#' @param x a \code{Ranges} object
#' @param shift shift how much should each range in x be shifted by? (default = 0L)
#' @param width width how long should the returned coverage score be?
#' This must be either a  positive integer, NA or NULL (default = NULL)
#' @param weight weight how much weight should be assigned to each range? Either
#' an integer or numeric vector or a column in x. (default = 1L)
#' @param method method see \link[IRanges]{coverage}.
#'
#' @return An expanded Ranges object with a score column corresponding to
#' the coverage value over that interval. Note that compute_coverage
#' drops metadata associated with the orginal ranges.
#' @seealso \link[IRanges]{coverage}
#' @importFrom IRanges coverage
#' @export
#' @rdname compute_coverage
#' @name compute_coverage.rd
setGeneric("compute_coverage",
           function(x, shift = 0L, width = NULL, weight = 1L,
                    method = c("auto", "sort", "hash")) {
             standardGeneric("compute_coverage")
             })

#' @describeIn compute_coverage compute the coverage of a GRanges object
setMethod("compute_coverage", "GRanges",
          function(x, shift = 0L, width = NULL, weight = 1L,
                   method = c("auto", "sort", "hash")) {
            GRanges(coverage(x, shift, width, weight, method))
          })

#' @describeIn compute_coverage compute the coverage of an IRanges object
setMethod("compute_coverage", "IRanges",
          function(x, shift = 0L, width = NULL, weight = 1L,
                   method = c("auto", "sort", "hash")) {
            cvg <- coverage(x, shift, width, weight, method)
            rng <- ranges(cvg)
            mcols(rng)[["score"]] <- runValue(cvg)
            rng
          })
