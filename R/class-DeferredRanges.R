setClassUnion("GenomicFile", 
              members = c("RsamtoolsFile", "RTLFile"))
# abstract class to represent operations performed on a file
setClass("FileOperator", 
         slots = c(input = "GenomicFile"),
         contains = "VIRTUAL")

setClass("BamFileOperator", 
         slots = c(
           param = "ScanBamParam",
           paired = "logical"
         ),
         contains = "FileOperator")

new_bam_op_list <- function(file, index = file, paired = FALSE) {
  if (is.null(index)) {
    index <- character()
  }
  new("BamOperationList",
             input = Rsamtools::BamFile(file, index = index),
             param = Rsamtools::ScanBamParam(),
             paired = paired)
}

# Deferred GRanges class

#' DeferredGenomicRanges Class
#'
#' @description  Enables deferred reading of BAM files by caching
#' after performing file operations.
#'
#' @param cache a GRanges object
#' @param operation an environment containing operations to perform
#' to generate the cache.
#' @param .data a DeferredGenomicRanges object
#'
#' @details
#' To see what is loaded in the cache use
#' `get_cache` and to see the file operations being performed
#' use `get_operation`. To check whether the cache is full
#' use `has_cache_contents`.
#' By default the `show` method for a DeferredGenomicRanges object
#' shows what's in the cache.
#'
#' @return a DeferredGenomicRanges object.
#' @seealso read_bam
setClass("DeferredGenomicRanges",
         slots = c(file_operation = "FileOperator"),
         contains = "DelegatingGenomicRanges")
