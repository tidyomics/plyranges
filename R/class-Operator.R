#' @importFrom Rsamtools BamFile ScanBamParam 
#' @export
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
load_genomic_file <- function(ops) UseMethod("load_genomic_file")

check_which_label <- function(ops) {
    length(bamWhich(ops@param)) > 0
}

paired_alignments <- function(ops) {
    alignments <- GenomicAlignments::readGAlignmentPairs(
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
    return(group_by(grng, "read_pair_id"))
}

unpaired_alignments <- function(ops) {
    galn_to_grng(
        GenomicAlignments::readGAlignments(file = ops@input,
                                       param = ops@param,
                                       with.which_label = check_which_label(ops))
    )
}

#' @importFrom Rsamtools bamFlag bamWhich ScanBamParam
#' @importFrom GenomicAlignments readGAlignments readGAlignmentPairs
load_genomic_file.BamFileOperator <- function(ops) {
    if (ops@paired) return(paired_alignments(ops))
    unpaired_alignments(ops)
}

#' @method select GRangesDeferred
#' @importFrom Rsamtools bamWhat<- bamTag<- scanBamWhat
#' @export
select.BamFileOperator <- function(.data, ..., .drop_ranges = FALSE) {
  dots <- quos(...)

  # populate bam params
  all_fields_tags <- unlist(lapply(dots, quo_name))

  tags_inx <- vapply(all_fields_tags,
                     FUN = function(x) { nchar(x) == 2L},
                     FUN.VALUE = logical(1))

  tags <- all_fields_tags[tags_inx]
  if (length(tags) > 0 ) bamTag(.data@param) <- tags

  fields <- all_fields_tags[!tags_inx]

  if (length(fields) > 0) {
    if (any(!(fields %in% Rsamtools::scanBamWhat()))) {
      stop("Invalid field identifier.", call. = FALSE)
    }
    bamWhat(.data@param) <- fields
  }
  return(.data)
}

# -- filter method for BamFileOperator
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

filter.BamFileOperator <- function(.data, ...) {
    dots <- quos(...)
    flags <- unlist(lapply(dots, quo_name))
    filters <- unlist(valid_flag_filters())
    flags <- intersect(names(filters), flags)
    if (length(flags) == 0L) {
        return(.data)
    } else {
        args <- as.list(rep(TRUE, length(flags)))
        names(args) <- unlist(lapply(flags, function(x) filters[[x]]))
        sam_flags <- do.call("scanBamFlag", args, envir = getNamespace("Rsamtools"))
    bamFlag(.data@ops$param) <- sam_flags
    }
}
