# Deferred GRanges class

#' GRangesDeferred  class
#'
#' @description  Enables deferred reading of BAM files
#'
#' @param cache a GRanges object
#' @param operation an environment containing operations to perform
#'
#' @export
#' @rdname ranges-deferred
setClass("GRangesDeferred",
         slots = c(cache = "GenomicRanges",
                   operation = "environment"),
         contains = "GRanges")

#' @rdname ranges-deferred
#' @export
GRangesDeferred <- function(cache = GRanges(), operation) {
  new("GRangesDeferred", cache = cache, operation = operation)
}

#' @export
setMethod("show", "GRangesDeferred",
          function(object) {
            show(get_cache(object))
          })
# -- accessors and checkers
#' @export
has_cache_contents <- function(.data) {
  length(.data@cache) > 0
}

#' @export
get_cache <- function(.data) { .data@cache }

#' @export
get_operation <- function(.data) { .data@operation }

