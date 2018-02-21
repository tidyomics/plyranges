#' @importFrom Rsamtools bamWhat<- bamTag<-
#' @export
select.GRangesDeferred <- function(.data, ..., .drop_ranges = FALSE) {
  dots <- quos(...)

  if (any(names(dots) != "")) {
    stop("Invalid input: ... must not be named.")
  }

  # if cached select columns
  if (has_cache_contents(.data)) {
    cache <- get_cache(.data)
    # no selection? - return .data
    if (length(dots) == 0) {
      return(cache)
    }
    return(select(cache, UQS(dots), .drop_ranges))
  }

  # else, populate bam params
  all_fields_tags <- unlist(lapply(dots, quo_name))

  tags_inx <- vapply(all_fields_tags,
                     FUN = function(x) { nchar(x) == 2L},
                     FUN.VALUE = logical(1))

  tags <- all_fields_tags[tags_inx]

  bamTag(.data@operation$param) <- tags

  fields <- all_fields_tags[!tags_inx]

  if (length(fields) > 0) {
    if (any(!(fields %in% Rsamtools::scanBamWhat()))) {
      stop("Invalid field identifier.", call. = FALSE)
    }
  }

  bamWhat(.data@operation$param) <- fields
  GRangesDeferred(load_alignments(.data),
                  operation = get_operation(.data))
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

#' @importFrom Rsamtools scanBamFlag bamFlag<-
#' @export
filter.GRangesDeferred <- function(.data, ...) {
  dots <- quos(...)

  # if cached normal filter
  if (has_cache_contents(.data)) {
    cache <- get_cache(.data)
    # no selection? - return .data
    if (length(dots) == 0) {
      return(cache)
    }
    return(filter(cache, UQS(dots)))
  }


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
    GRangesDeferred(load_alignments(.data),
                    operation = get_operation(.data))
  }
}



#' @importFrom Rsamtools bamWhich<-
#' @export
filter_by_overlaps.GRangesDeferred <- function(x, y, maxgap = -1L, minoverlap = 0L) {

  stopifnot(is(y, "GenomicRanges"))

  if (has_cache_contents(x)) {
    cache <- get_cache(x)
    return(filter_by_overlaps(cache, y, maxgap = -1L, minoverlap = 0L))
  }

  bamWhich(x@operation$param) <- y

  GRangesDeferred(cache = load_alignments(x),
                  operation = get_operation(x))
}

#' @export
summarise.GRangesDeferred <- function(.data, ...) {
  dots <- quos(...)
  if (has_cache_contents(.data)) {
    cache <- get_cache(.data)
  } else {
    cache <- load_alignments(.data)
  }
  return(summarize(cache, UQS(dots)))
}

#' @export
mutate.GRangesDeferred <- function(.data, ...) {
  dots <- quos(...)
  if (has_cache_contents(.dat)) {
    cache <- get_cache(.data)
  } else {
    cache <- load_alignments(.data)
  }
  return(mutate(cache, UQS(dots)))
}

#' @export
group_by.GRangesDeferred <- function(.data, ...) {
  capture_groups <- quos(...)
  groups <- lapply(capture_groups, function(x) quo_name(x))
  groups <- syms(groups)

  if (has_cache_contents(.data)) {
    cache <- get_cache(.data)
  } else {
    cache <- load_alignments(.data)
  }
  new("GRangesGrouped", cache, groups = groups)
}

