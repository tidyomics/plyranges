#' Read a BAM file
#'
#' @param file A connection or path to a BAM file
#' @param index The path to the BAM index file
#' @param paired Whether to treat alignments as paired end (TRUE)
#' or single end (FALSE). Default is FALSE.
#'
#' @details Reading a BAM file is deferred until an action
#' such as using `summarise()` or `mutate()` occurs. If paired is set to
#' TRUE, when alignments are loaded, the GRanges has two additional
#' columns called read_pair_id and read_pair_group corresponding
#' to paired reads and is grouped by the read_pair_group.
#'
#' Certain verbs have different behaviour, after using `read_bam()`.
#' 
#' For `select()` valid columns are the fields available in the
#' BAM file. Valid entries are qname (QNAME), flag (FLAG),
#' rname (RNAME), strand, pos (POS), qwidth (width of query),
#' mapq (MAPQ), cigar (CIGAR), mrnm (RNEXT), mpos (PNEXT), isize
#' (TLEN), seq (SEQ), and qual (QUAL). Any two character tags in the BAM file 
#' are also valid.
#'
#' For `filter()` the following fields are valid, to select the FALSE option
#' place `!` in front of the field:
#' 
#' * `is_paired` Select either unpaired (FALSE) or paired (TRUE) reads.
#' * `is_proper_pair` Select either improperly paired (FALSE) or properly
#' paired (TRUE) reads. This is dependent on the alignment software used.
#' * `is_unmapped_query``	Select unmapped (TRUE) or mapped (FALSE) reads.
#' * `has_unmapped_mate` Select reads with mapped (FALSE) or unmapped (TRUE) mates.
#' * `is_minus_strand` 	Select reads aligned to plus (FALSE) or minus (TRUE) strand.
#' * `is_mate_minus_strand`	Select reads where mate is aligned to plus (FALSE) or
#' minus (TRUE) strand.
#' * `is_first_mate_read`	Select reads if they are the first mate (TRUE) or
#' not (FALSE).
#' * `is_second_mate_read` Select reads if they are the second mate (TRUE) or
#' not (FALSE).
#' * `is_secondary_alignment` Select reads if their alignment status is
#' secondary (TRUE) or not (FALSE). This might be relevant if there are
#' multimapping reads.
#' * `is_not_passing_quality_controls` Select reads that either pass
#' quality controls (FALSE) or that do not (TRUE).
#' * `is_duplicate`	Select reads that are unduplicated (FALSE) or
#' duplicated (TRUE). This may represent reads that are PCR or
#' optical duplicates.
#'
#' @importFrom GenomicAlignments readGAlignments readGAlignmentPairs
#' @importFrom Rsamtools BamFile ScanBamParam
#' @export
#' @examples
#' if (require(pasillaBamSubset)) {
#'    bamfile <- untreated1_chr4()
#'    # nothing is read until an action has been performed
#'    print(read_bam(bamfile))
#'    # define a region of interest
#'    roi <- data.frame(seqnames = "chr4", start = 5e5, end = 7e5) %>%
#'             as_granges()
#'    rng <- read_bam(bamfile) %>% 
#'             select(mapq) %>%
#'             filter_by_overlaps(roi)
#' }
#'
#' @return A DeferredGenomicRanges object
#' @rdname io-bam-read
#' @seealso [Rsamtools::BamFile()],[GenomicAlignments::readGAlignments()]
read_bam <- function(file, index = file, paired = FALSE) {
  ops <- new_bam_ops(file, index, paired)
  new_DeferredGenomicRanges(GRanges(), ops)
}