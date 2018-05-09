#' Group a GRanges object by introns or gaps
#' 
#' @param x a GenomicRanges object with a cigar string column
#' 
#' @details Creates a grouped Ranges object from a cigar string
#' column, for `chop_by_introns()` will check for the presence of
#' "N" in the cigar string and create a new column called
#' `intron` where TRUE indicates the alignment has a skipped
#' region from the reference. For `chop_by_gaps()` will check
#' for the presence of "N" or "D" in the cigar string and
#' create a new column called "gaps" where TRUE indicates
#' the alignment has a deletion from the reference or has an intron.  
#' 
#' @importFrom S4Vectors Rle
#' @rdname ranges-chop
#' @export
chop_by_introns <- function(x) UseMethod("chop_by_introns")

#' @export
chop_by_introns.GenomicRanges <- function(x) {
  # check for cigar column
  stopifnot(any(names(mcols(x)) %in% "cigar"))
  # add grouping var
  mcols(x)$intron <- Rle(grepl("N", mcols(x)$cigar, fixed = TRUE))
  group_by(x, "intron")
}

chop_by_introns.DeferredGenomicRanges <- function(x) {
  chop_by_introns(load_delegate(x))
}

#' @rdname ranges-chop
#' @export
chop_by_gaps <- function(x) UseMethod("chop_by_gaps")

chop_by_gaps.GenomicRanges <- function(x) {
  # check for cigar column
  stopifnot(any(names(mcols(x)) %in% "cigar"))
  # add grouping var
  mcols(x)$gaps <- Rle(grepl("N|D", mcols(x)$cigar))
  group_by(x, "gaps")
}

chop_by_gaps.DeferredGenomicRanges <- function(x) {
  chop_by_gaps(load_delegate(x))
}