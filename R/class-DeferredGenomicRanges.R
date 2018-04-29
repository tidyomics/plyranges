
#' File Operator class
#' @description A virtual class that represents actions to be performed
#' on common data formats.
#' @rdname ranges-file-operator
#' @export
setClass("FileOperator", contains = "VIRTUAL")

# A class for performing deferred reading of a file
#' DeferredGenomicRanges Class
#'
#' @description  Enables deferred reading of files (currently only BAM files) by 
#' caching results after a plyranges verb is called.
#'
#' @slot delegate a GenomicRanges object to be cached
#' @slot ops  A FileOperator object
#'
#' @details
#' To see what is loaded in the cache use
#' `get_cache` and to see the file operations being performed
#' use `get_operation`. To check whether the cache is full
#' use `has_cache_contents`.
#' By default the `show` method for a GRangesDeferred object
#' shows what's in the cache.
#'
#' @return a GRangesDeferred object.
#' @seealso read_bam
#' @export
#' @rdname ranges-deferred
setClass("DeferredGenomicRanges", 
          slots = c(ops = "FileOperator"),
          contains = "DelegatingGenomicRanges"
)

new_DeferredGenomicRanges <- function(delegate, ops) {
    new("DeferredGenomicRanges", delegate = delegate, ops = ops)
}

is_empty_delegate <- function(.data) length(.data@delegate) == 0L

# --- dplyr verbs for deferred classes
#' #' @method select DeferredGenomicRanges
#' #' @importFrom Rsamtools bamWhat<- bamTag<-
#' #' @export
#' select.GRangesDeferred <- function(.data, ..., .drop_ranges = FALSE) {
#'   dots <- quos(...)
#' 
#'   if (any(names(dots) != "")) {
#'     stop("Invalid input: ... must not be named.")
#'   }
#'   # if no cache update ops
#'   if (is_empty_delegate(.data)) {
#'       ops <- select(.data@ops, UQS(dots), .drop_ranges)
#'       new_DeferredGenomicRanges(delegate = load_genomic_file(ops), ops = ops)
#'   } else {
#'     return(select(.data, UQS(dots), .drop_ranges))  
#'   }
#' }