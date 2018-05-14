#' @method select BamFileOperator
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

#' @importFrom Rsamtools bamMapqFilter<-
filter.BamFileOperator <- function(.data, ...) {
  dots <- quos(...)
  filters <- unlist(valid_flag_filters())
  
  flags <- unlist(lapply(dots, quo_name))
  
  # mapq filter 
  mapq_filter <- grepl("mapq", flags)
  mapq_args <- flags[mapq_filter]
  if (length(mapq_args) > 0) {
    mapq_vals <- as.integer(sub(".*>", "", mapq_args))
    if (any(is.na(mapq_vals))) {
      message("mapq filter only accepts greater than statements")
    }
    mapq_vals <- mapq_vals[!is.na(mapq_vals)]
    if (length(mapq_vals) > 1) {
      message("mapq filter only accepts singular mapq filters; taking first
              argument")
    } 
    bamMapqFilter(.data@param) <- mapq_vals
  }
  
  # flag filters
  args <- !grepl("^!", flags[!mapq_filter])
  names(args) <- gsub("^!", "", flags[!mapq_filter])
  
  valid_flags <- intersect(names(filters), names(args))
  
  if (length(valid_flags) == 0L & length(mapq_args) == 0L) {
    stop("no valid flags found in filter", call. = FALSE)
  }
  
  if (any(!(names(args) %in% valid_flags))) {
    invalid_flag <- names(args)[!(names(args) %in% valid_flags)]
    stop(paste(paste(invalid_flag, collapse = ","), "are not valid flags"),
         call. = FALSE)
  }
  names(args) <- filters[valid_flags]
  args <- as.list(args)
  sam_flags <- do.call("scanBamFlag", 
                       args, 
                       envir = getNamespace("Rsamtools"))
  bamFlag(.data@param) <- sam_flags
  .data
}

#' @importFrom Rsamtools bamWhich<-
filter_by_overlaps.BamFileOperator <- function(x, y, maxgap = -1L, minoverlap = 0L) {
  stopifnot(is(y, "GenomicRanges"))
  bamWhich(x@param) <- y
  x
}
