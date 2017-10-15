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

#' @importFrom Rsamtools scanBamFlag
bam_flag_reads <- function(is_paired = NULL,
                           is_proper_pair = NULL,
                           is_unmapped_query = NULL,
                           has_unmapped_mate = NULL,
                           is_minus_strand = NULL,
                           is_mate_minus_strand = NULL,
                           is_first_mate_read = NULL,
                           is_second_mate_read = NULL,
                           is_not_primary_read = NULL,
                           is_secondary_alignment = NULL,
                           is_not_passing_quality_control = NULL,
                           is_duplicate = NULL) {

  cl <- Filter(function(.) !is.symbol(.), as.list(match.call()))
  if (length(cl) > 0) {
    cl <- Map(check_lgl, cl)
    names(cl) <- to_camel(names(cl))
  }


  if (requireNamespace("Rsamtools")) {
    scanBamFlag <- Rsamtools::scanBamFlag
  } else {
    stop("Package Rsamtools must be installed to use this function.",
         call. = FALSE)
  }
  do.call("scanBamFlag", cl)

}

#' @importFrom Rsamtools scanBamParam
bam_select_reads <- function(fields = NULL, tags = NULL, overlap_ranges = NULL) {

  if (is.null(fields)) {
    fields <- character(0)
  } else {
    if (!(fields %in% Rsamtools::scanBamWhat())) {
      stop("Invalid field identifier.", call. = FALSE)
    }
  }
  if (is.null(tags)) {
    tags <- character(0)
  }

  if (!is.null(overlap_ranges)) {
    stopifnot(is(overlap_ranges, "GRanges"))
  } else {
    overlap_ranges <- GRanges()
  }

  if (requireNamespace("Rsamtools")) {
    Rsamtools::ScanBamParam(what = fields,
                            tag = tags,
                            which = overlap_ranges)
  } else {
    stop("Package: Rsamtools must be installed to use this function.",
         call. = FALSE)
  }
}


read_bam <- function(file, index = file, flag_reads = bam_flag_reads(),
                     select_reads = bam_select_reads(),
                     paired = TRUE) {

  if (paired) {
    io_bam <- GenomicAlignments::readGAlignmentPairs
  } else {
    io_bam <- GenomicAlignments::readGAlignments
  }
  param <- select_reads()
  Rsamtools::bamFlag(param) <- flag_reads

  alignments <- io_bam(file, index, param = param)
  as(alignments, "GRanges")
}
