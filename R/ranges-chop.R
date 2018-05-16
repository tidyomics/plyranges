expand_rng_by_cigar <- function(x, type) {
  # check for cigar column
  stopifnot(any(names(mcols(x)) %in% "cigar"))
  if (type == "gaps") {
    drop.D.ranges <- TRUE
  } else {
    drop.D.ranges <- FALSE
  }
  # extract alignment ranges, returns an IRangesList
  rng <- extractAlignmentRangesOnReference(mcols(x)$cigar, 
                                           pos=start(x), 
                                           drop.D.ranges = drop.D.ranges)
  
  n <- elementNROWS(rng)
  seqnames <- BiocGenerics::rep.int(seqnames(x), n)
  strand <- BiocGenerics::rep.int(strand(x), n)
  grng <- GRanges(seqnames = seqnames, 
                  strand = strand, 
                  ranges = rng@unlistData)
  
  seqinfo(grng) <- seqinfo(x)
  grp <- BiocGenerics::rep.int(seq_along(x), n)
  mcols(grng) <- mcols(x)[grp,]
  mcols(grng)[[type]] <- grp
  group_by(grng, UQ(type))
}

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
#' @importFrom BiocGenerics rep.int
#' @importFrom S4Vectors Rle elementNROWS
#' @importFrom GenomicAlignments extractAlignmentRangesOnReference
#' @rdname ranges-chop
#' @export
chop_by_introns <- function(x) UseMethod("chop_by_introns")

#' @export
chop_by_introns.GenomicRanges <- function(x) {
  expand_rng_by_cigar(x, "introns")
}
#' @export
chop_by_introns.DeferredGenomicRanges <- function(x) {
  chop_by_introns(load_delegate(x))
}

#' @rdname ranges-chop
#' @export
chop_by_gaps <- function(x) UseMethod("chop_by_gaps")

#' @export
chop_by_gaps.GenomicRanges <- function(x) {
  expand_rng_by_cigar(x, "gaps")
}

#' @export
chop_by_gaps.DeferredGenomicRanges <- function(x) {
  chop_by_gaps(load_delegate(x))
}