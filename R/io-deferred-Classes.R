# Deferred GRanges class

#' GRangesDeferred Class
#'
#' @description  Enables deferred reading of BAM files by caching
#' after performing file operations.
#'
#' @param cache a GRanges object
#' @param operation an environment containing operations to perform
#' to generate the cache.
#'
#' @details
#' To see what is loaded in the cache use
#' \code{get_cache} and to see the file operations being performed
#' use \code{get_operation}. To check whether the cache is full
#' use \code{has_cache_contents}.
#' By default the \code{show} method for a GRangesDeferred object
#' shows what's in the cache.
#'
#' @seealso read_bam
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

#' @rdname ranges-deferred
setMethod("show", "GRangesDeferred",
          function(object) {
            show(get_cache(object))
          })
# -- accessors and checkers

#' @rdname ranges-deferred
#' @export
has_cache_contents <- function(.data) {
  length(.data@cache) > 0
}

#' @rdname ranges-deferred
#' @export
get_cache <- function(.data) { .data@cache }

#' @rdname ranges-deferred
#' @export
get_operation <- function(.data) { .data@operation }

