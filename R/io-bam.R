# io-bam

#' Constructor for deferred reading of BAM files
#'
#' @param file The path to a BAM file
#' @param index The path to the BAM index file
#' @importFrom Rsamtools BamFile ScanBamParam
#' @export
bam_file <- function(file, index = file) {

  env <- rlang::new_environment(
    list(input = Rsamtools::BamFile(file, index = index),
         param = Rsamtools::ScanBamParam())
    )
  GRangesDeferred(env)
}


#' @importFrom Rsamtools bamTag<-
#' @export
select_tags <- function(bam, ...) {
  stopifnot(is(bam, "GRangesDeferred"))

  tags <- quos(...)
  if (any(names(tags) != "")) {
    stop("Invalid input: ... must not be named.")
  }


  if (length(tags) == 0) {
    tags <- character(0)
  } else {
    # check two character long input
    valid_tag <- all(vapply(tags, nchar, FUN.VALUE = integer(1)) == 2L)
    if (!valid_tag) {
      stop("Input tags must have length 2.")
    }
  }

  bamTag(bam@operation$param) <- tags
  return(bam)
}

#' Select fields and tags to load from a BAM
#'
#' @param bam a BamFile
#' @param ... Values indicating which fields from
#' the BAM file to load. Valid entries are qname (QNAME), flag (FLAG),
#' rname (RNAME), strand, pos (POS), qwidth (width of query),
#' mapq (MAPQ), cigar (CIGAR), mrnm (RNEXT), mpos (PNEXT), isize
#' (TLEN), seq (SEQ), qual (QUAL).
#' @seealso \code{\link[Rsamtools]{ScanBamParam}}
#' @importFrom Rsamtools bamWhat<-
#' @export
select_fields <- function(bam, ...) {

  stopifnot(is(bam, "GRangesDeferred"))

  fields <- quos(...)
  if (any(names(fields) != "")) {
    stop("Invalid input: ... must not be named.")
  }

  fields <- unlist(lapply(fields, quo_name))
  if (length(fields) == 0L) {
    fields <- character(0)
  } else {
    if (any(!(fields %in% Rsamtools::scanBamWhat()))) {
      stop("Invalid field identifier.", call. = FALSE)
    }
    bamWhat(bam@operation$param) <- fields
  }
  bam
}

valid_flag_filters <- function() {
  list(
    is_paired = "isPaired",
    is_proper_pair = "isProperPair",
    is_unmapped_query = "isUnmappedQuery",
    has_unmapped_mate = "hasUnmappedMate",
    is_minus_strand = "isMinusStrand",
    is_mate_minus_strand = "isMateMinusStrand",
    is_first_mate_read = "isFirstMateRead",
    is_second_mate_read = "isSecondMateRead",
    is_secondary_alignment = "isSecondaryAlignment",
    is_not_passing_quality_controls = "isNotPassingQualityControls",
    is_duplicate = "isDuplicate"
  )

}
#' Flag reads for filtering prior to loading a BAM file.
#'
#' @param .data a DeferredGRanges object
#' @param ... filter to apply to the input file
#'
#' @details The following are valid filters for a `BamFile`:
#' \describe{
#'  \item{is_paired}	Select either unpaired (FALSE) or paired (TRUE) reads.
#'  \item{is_proper_pair} Select either improperly paired (FALSE) or properly paired (TRUE) reads. This is dependent on the alignment software used.
#'  \item{is_unmapped_query}	Select unmapped (TRUE) or mapped (FALSE) reads.
#'  \item{has_unmapped_mate} Select reads with mapped (FALSE) or unmapped (TRUE)  mates.
#'  \item{is_minus_strand}	Select reads aligned to plus (FALSE) or minus (TRUE) strand.
#'  \item{is_mate_minus_strand}	Select reads where mate is aligned to plus (FALSE) or minus (TRUE) strand.
#'  \item{is_first_mate_read}	Select reads if they are the first mate (TRUE) or not (FALSE).
#'  \item{is_second_mate_read}	Select reads if they are the second mate (TRUE) or not (FALSE).
#'  \item{is_secondary_alignment Select reads if their alignment status is secondary (TRUE) or not (FALSE). This might be relevant if there are multimapping reads.
#'  \item{is_not_passing_quality_controls} Select reads that either pass quality controls (FALSE) or that do not (TRUE).
#'  \item{is_duplicate}	Select reads that are unduplicated (FALSE) or duplicated (TRUE). This may represent reads that are PCR or optical duplicates.
#' }
#'
#' @importFrom Rsamtools scanBamFlag bamFlag<-
filter.GRangesDeferred <- function(.data, ...) {

  dots <- quos(...)
  flags <- unlist(lapply(dots, quo_name))
  filters <- unlist(valid_flag_filters())
  flags <- intersect(names(filters), flags)

  if (length(flags) == 0L) {
    return(.data)
  } else {
    args <- as.list(rep(TRUE, length(flags)))
    names(args) <- unlist(lapply(flags, function(x) filters[[x]]))
    sam_flags <- do.call("scanBamFlag",
                         args,
                         envir = getNamespace("Rsamtools"))
    bamFlag(.data@operation$param) <- sam_flags
    return(.data)
  }
}

#' @importFrom Rsamtools bamWhich<-
filter_by_overlaps.GRangesDeferred <- function(x, y, maxgap = NA, minoverlap = NA) {
  stopifnot(is(y, "GenomicRanges"))
  bamWhich(x@operation$param) <- y
  x
}

galn_to_grng <- function(alignments) {
  grng <- granges(alignments, use.mcols = TRUE)
  mcols(grng)[["cigar"]] <- GenomicAlignments::cigar(alignments)
  mcols(grng)[["qwidth"]] <- GenomicAlignments::qwidth(alignments)
  mcols(grng)[["njunc"]] <- GenomicAlignments::njunc(alignments)
  grng
}
#' Read a BAM file
#'
#' @param bam A GRangesDeferred object (constructed with \code{bam_file}.
#' @param paired Whether to treat the BAM file as paired end (TRUE)
#' or single end (FALSE).
#' @details By default for paired = FALSE, read_bam will select reads
#' that are mapped, paired and have mapped mates. For paired = FALSE,
#' read_bam will select reads that are mapped. This function is wrapper to
#' the \code{readGAlignment} functions in \pkg{GenomicAlignments}. To perform
#' any filtering or selection prior to filtering use the \code{select_tags},
#' \code{select_fields} or \select{filter} methods provided.
#'
#' @seealso \code{\link[GenomicAlignments]{readGAlignments}}
#' @importFrom Rsamtools bamFlag bamWhich ScanBamParam
#' @importFrom GenomicAlignments readGAlignments readGAlignmentPairs
#' @export
read_bam <- function(bam, paired = FALSE) {

  if (is(bam, "GRangesDeferred")) {
    if (length(bamWhich(bam@operation$param)) > 0) {
      with.which_label <- TRUE
    } else {
      with.which_label <- FALSE
    }
  } else {
    stop("Unable to read BAM.", call. = FALSE)
  }

  file <- bam@operation$input
  param <- bam@operation$param

  if (paired) {
    alignments <- GenomicAlignments::readGAlignmentPairs(
      file,
      param = param,
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
    return(group_by(grng, "read_pair_group"))

  } else {
    alignments <- GenomicAlignments::readGAlignments(
      file,
      param = param,
      with.which_label = with.which_label
    )
    return(galn_to_grng(alignments))
  }
}
