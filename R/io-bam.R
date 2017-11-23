# io-bam
check_lgl <- function(x) {
  stopifnot(is.logical(x) || is.na(x))
  x
}

to_camel <- function(x) {
  capit <- function(x) paste0(toupper(substring(x, 1, 1)),
                              substring(x, 2, nchar(x)))

  sapply(strsplit(x, "\\_"),
         function(x) paste0(x[1], paste(capit(x[-1]), collapse = "")))
}

#' Flag reads for filtering prior to loading a BAM file.
#'
#' @param is_paired	Select either unpaired (FALSE) or paired (TRUE) reads.
#' @param is_proper_pair Select either improperly paired (FALSE) or properly paired (TRUE) reads. This is dependent on the alignment software used.
#' @param is_unmapped_query	Select unmapped (TRUE) or mapped (FALSE) reads.
#' @param has_unmapped_mate Select reads with mapped (FALSE) or unmapped (TRUE)  mates.
#' @param is_minus_strand	Select reads aligned to plus (FALSE) or minus (TRUE) strand.
#' @param is_mate_minus_strand	Select reads where mate is aligned to plus (FALSE) or minus (TRUE) strand.
#' @param is_first_mate_read	Select reads if they are the first mate (TRUE) or not (FALSE).
#' @param is_second_mate_read	Select reads if they are the second mate (TRUE) or not (FALSE).
#' @param is_secondary_alignment Select reads if their alignment status is secondary (TRUE) or not (FALSE). This might be relevant if there are multimapping reads.
#' @param is_not_passing_quality_control Select reads that either pass quality controls (FALSE) or that do not (TRUE).
#' @param is_duplicate	Select reads that are unduplicated (FALSE) or duplicated (TRUE). This may represent reads that are PCR or optical duplicates.
#'
#' @details By default all arguments are set to NULL which means any read can
#' be selected and their status according to the above parameters are ignored.
#'
#' @seealso \code{\link[Rsamtools]{ScanBamParam}}
#' @importFrom Rsamtools scanBamFlag
#' @export
bam_flag_reads <- function(is_paired = NULL,
                           is_proper_pair = NULL,
                           is_unmapped_query = NULL,
                           has_unmapped_mate = NULL,
                           is_minus_strand = NULL,
                           is_mate_minus_strand = NULL,
                           is_first_mate_read = NULL,
                           is_second_mate_read = NULL,
                           is_secondary_alignment = NULL,
                           is_not_passing_quality_control = NULL,
                           is_duplicate = NULL) {

  cl <- Filter(function(.) !is.symbol(.), as.list(match.call()))
  if (length(cl) > 0) {
    cl <- Map(check_lgl, cl)
    names(cl) <- to_camel(names(cl))
  }

  do.call("scanBamFlag", cl, envir = getNamespace("Rsamtools"))

}

#' Select fields and tags to load from a BAM
#'
#' @param fields a character vector indicating which fields from
#' the BAM file to load. Valid entries are "qname" (QNAME), "flag" (FLAG),
#' "rname" (RNAME), "strand", "pos" (POS), "qwidth" (width of query),
#' "mapq" (MAPQ), "cigar" (CIGAR), "mrnm" (RNEXT), "mpos" (PNEXT), "isize"
#' (TLEN), "seq" (SEQ), "qual" (QUAL).
#' @param tags a character vector of tags (optional fields associated
#' with each read), these are two letter codes.
#' @seealso \code{\link[Rsamtools]{ScanBamParam}}
#' @importFrom Rsamtools ScanBamParam
#' @export
bam_select <- function(fields = NULL, tags = NULL) {

  if (is.null(fields)) {
    fields <- character(0)
  } else {
    if (any(!(fields %in% Rsamtools::scanBamWhat()))) {
      stop("Invalid field identifier.", call. = FALSE)
    }
  }
  if (is.null(tags)) {
    tags <- character(0)
  } else {
    # check two character long input
    stopifnot(all(vapply(tags, nchar, FUN.VALUE = integer(1)) == 2L))
  }

  Rsamtools::ScanBamParam(what = fields,
                          tag = tags)
}


#' Read a BAM file
#'
#' @param file The path to a BAM file.
#' @param index  The path to the index file of the BAM. Must be given without
#' the '.bai' extension.
#' @param paired Whether to treat the BAM file as paired end (TRUE) or single end (FALSE.
#' @param flag_reads The reads that will be returned by selecting combinations in
#'\code{bam_flag_reads}.
#' @param select_reads The fields and tags in the BAM file that will be returned
#' as metadata columns in the returned GRanges object. Determined by \code{bam_select}.
#' @param overlap_ranges Only return records that overlap the given Ranges object. The
#' returning GRanges object will be annotated with a column called which_label to
#' identify the range in overlap_ranges that the read overalps.
#'
#' @details By default for paired = TRUE, read_bam will select reads
#' that are mapped, paired and have mapped mates. For paired = FALSE,
#' read_bam will select reads that are mapped. This function is wrapper to
#' the \code{readGAlignment} functions in \pkg{GenomicAlignemnts}.
#' @seealso \code{\link[GenomicAlignments]{readGAlignments}}
#' @importFrom Rsamtools bamFlag bamWhich ScanBamParam
#' @importFrom GenomicAlignments readGAlignments readGAlignmentPairs
#' @export
read_bam <- function(file, index = file,
                     paired = TRUE,
                     flag_reads = bam_flag_reads(),
                     select_reads = bam_select(),
                     overlap_ranges = NULL
                     ) {

  if (paired) {
    io_bam <- GenomicAlignments::readGAlignmentsList
  } else {
    io_bam <- GenomicAlignments::readGAlignments
  }

  if (!is.null(overlap_ranges)) {
    stopifnot(is(overlap_ranges, "GRanges"))
    with.which_label <- TRUE
  } else {
    overlap_ranges <- GRanges()
    with.which_label <- FALSE
  }

  param <- select_reads
  Rsamtools::bamFlag(param) <- flag_reads
  Rsamtools::bamWhich(param) <- overlap_ranges

  alignments <- io_bam(file, index,
                       param = param,
                       with.which_label = with.which_label)
  alignments
}
