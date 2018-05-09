#' @include class-Operator.R

#' @title DeferredGenomiRanges objects
#'
#' @description  Enables deferred reading of files (currently only BAM files) by 
#' caching results after a plyranges verb is called.
#'
#' @slot delegate a GenomicRanges object to be cached
#' @slot ops  A FileOperator object
#' 
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

load_delegate <-   function(.data) {
  if (is_empty_delegate(.data)) {
    load_genomic_file(.data@ops)
  } else {
    .data@delegate
  }
}