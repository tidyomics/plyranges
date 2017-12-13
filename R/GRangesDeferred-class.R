#' Deferred GRanges class
#'
#'
setClass("GRangesDeferred",
         slots = c(cache = "GenomicRanges",
                   operation = "environment"),
         contains = "GRanges")


#' Constructor for GRangesDeferred class
#'
#' @param operation an environment
#' @export
GRangesDeferred <- function(operation) {
  if (typeof(operation) != "environment") {
    stop("input must be an environment")
  }
  cache <- GRanges()
  new("GRangesDeferred", cache = cache, operation = operation)
}
