# abstract class to represent operations 
# performed over a File

#' An abstract class to represent operations performed over a file
#' 
#' @details This class is used internally by DeferredGenomicRanges objects.
#' Currently, this class is only implemented for bam files (as a 
#' BamFileOperator) but will eventually be extended to the other avaialable
#' readers.
#' 
#' @rdname ranges-class
#' @export
setClass("FileOperator", contains = "VIRTUAL")

#' @importFrom Rsamtools BamFile ScanBamParam 
#' @export
#' @rdname ranges-class
setClass("BamFileOperator", 
         slots = c(
           input = "BamFile",
           param = "ScanBamParam",
           paired = "logical"
         ),
         contains = "FileOperator"
)

new_bam_ops <- function(file, index = file, paired = FALSE) {
  if (is.null(index)) {
    index <- character()
  }
  new("BamFileOperator", 
      input = Rsamtools::BamFile(file, index = index),
      param = Rsamtools::ScanBamParam(),
      paired = paired
  )
}

# Generic File loader
#' @importFrom Rsamtools bamFlag bamWhich ScanBamParam
#' @importFrom GenomicAlignments readGAlignments readGAlignmentPairs
#' @importFrom S4Vectors first second
load_genomic_file <- function(ops) UseMethod("load_genomic_file")


load_genomic_file.BamFileOperator <- function(ops) {
  if (ops@paired) return(paired_alignments(ops))
  unpaired_alignments(ops)
}

check_which_label <- function(ops) {
  length(bamWhich(ops@param)) > 0
}

galn_to_grng <- function(alignments) {
  grng <- granges(alignments, use.mcols = TRUE)
  mcols(grng)[["cigar"]] <- GenomicAlignments::cigar(alignments)
  mcols(grng)[["qwidth"]] <- GenomicAlignments::qwidth(alignments)
  mcols(grng)[["njunc"]] <- GenomicAlignments::njunc(alignments)
  seqinfo(grng) <- seqinfo(alignments)
  grng
}

paired_alignments <- function(ops) {
  alignments <- readGAlignmentPairs(
    ops@input,
    param = ops@param,
    with.which_label = check_which_label(ops)
  )
  r1 <- first(alignments)
  r2 <- second(alignments)
  r1_grng <- galn_to_grng(r1)
  mcols(r1_grng)[["read_pair_id"]] <- seq_along(r1_grng)
  mcols(r1_grng)[["read_pair_group"]] <- 1L
  r2_grng <- galn_to_grng(r2)
  mcols(r2_grng)[["read_pair_id"]] <- seq_along(r1_grng)
  mcols(r2_grng)[["read_pair_group"]] <- 2L
  grng <- c(r1_grng, r2_grng)
  grng <- grng[order(grng$read_pair_id)]
  return(group_by(grng, !!rlang::sym("read_pair_id")))
}

unpaired_alignments <- function(ops) {
  galn_to_grng(
    suppressWarnings(
      readGAlignments(
        file = ops@input,
        param = ops@param,
        with.which_label = check_which_label(ops)
      )
    )
  )
}