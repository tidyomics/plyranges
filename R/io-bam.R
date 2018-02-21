#' Read a BAM file
#'
#' @param file A connection or path to a BAM file
#' @param index The path to the BAM index file
#' @param paired Whether to treat the BAM file as paired end (TRUE)
#' or single end (FALSE).
#'
#' @details Reading a BAM file is deferred until an action
#' such as selecting tags or filtering by overlapping ranges
#' or filtering by tags is performed on it. If paired is set to
#' TRUE, when alignments are load, the GRanges has two additional
#' columns called read_pair_id and read_pair_group corresponding
#' to paired reads (by default paired = TRUE will only select
#' reads that are proper pairs).
#'
#' For \code{select} valid columns are the either the fields of the
#' BAM file. Valid entries are qname (QNAME), flag (FLAG),
#' rname (RNAME), strand, pos (POS), qwidth (width of query),
#' mapq (MAPQ), cigar (CIGAR), mrnm (RNEXT), mpos (PNEXT), isize
#' (TLEN), seq (SEQ), and qual (QUAL). Any two character
#' tags in the BAM file are also valid.
#'
#' For \code{filter} the following fields are valid
#' is_paired Select either unpaired (FALSE) or paired (TRUE) reads.
#' is_proper_pair Select either improperly paired (FALSE) or properly
#' paired (TRUE) reads. This is dependent on the alignment software used.
#' is_unmapped_query	Select unmapped (TRUE) or mapped (FALSE) reads.
#' has_unmapped_mate Select reads with mapped (FALSE) or unmapped (TRUE) mates.
#' is_minus_strand 	Select reads aligned to plus (FALSE) or minus (TRUE) strand.
#' is_mate_minus_strand	Select reads where mate is aligned to plus (FALSE) or
#' minus (TRUE) strand.
#' is_first_mate_read	Select reads if they are the first mate (TRUE) or
#' not (FALSE).
#' is_second_mate_read Select reads if they are the second mate (TRUE) or
#' not (FALSE).
#' is_secondary_alignment Select reads if their alignment status is
#' secondary (TRUE) or not (FALSE). This might be relevant if there are
#' multimapping reads.
#' is_not_passing_quality_controls Select reads that either pass
#' quality controls (FALSE) or that do not (TRUE).
#' is_duplicate	Select reads that are unduplicated (FALSE) or
#' duplicated (TRUE). This may represent reads that are PCR or
#' optical duplicates.
#'
#'
#' @importFrom GenomicAlignments readGAlignments readGAlignmentPairs
#' @importFrom Rsamtools BamFile ScanBamParam
#' @importFrom rlang new_environment
#' @export
#' @rdname io-bam-read
read_bam <- function(file, index = file, paired = FALSE) {
  env <- rlang::new_environment(
    list(input = Rsamtools::BamFile(file, index = index),
         param = Rsamtools::ScanBamParam(),
         paired = paired)
    )
  GRangesDeferred(operation = env)
}


#' @importFrom Rsamtools bamFlag bamWhich ScanBamParam
#' @importFrom GenomicAlignments readGAlignments readGAlignmentPairs
load_alignments <- function(.data) {

  if (length(bamWhich(.data@operation$param)) > 0) {
    with.which_label <- TRUE
  } else {
    with.which_label <- FALSE
  }

  if (.data@operation$paired) {
    alignments <- GenomicAlignments::readGAlignmentPairs(
      .data@operation$input,
      param = .data@operation$param,
      with.which_label = with.which_label
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
    return(grng)
  }
  alignments <-  GenomicAlignments::readGAlignments(file = .data@operation$input,
                      param = .data@operation$param,
                      with.which_label = with.which_label)
  galn_to_grng(alignments)
}


galn_to_grng <- function(alignments) {
  grng <- granges(alignments, use.mcols = TRUE)
  mcols(grng)[["cigar"]] <- GenomicAlignments::cigar(alignments)
  mcols(grng)[["qwidth"]] <- GenomicAlignments::qwidth(alignments)
  mcols(grng)[["njunc"]] <- GenomicAlignments::njunc(alignments)
  grng
}

