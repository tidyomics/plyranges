# Representing seqinfo information as a Ranges

#' Construct annotation information for Ranges
#' @param seqnames A character vector containing the name of sequences.
#' @param seqlengths An optional integer vector containg the lengths of sequences.
#' @param is_circular An optional logical vector indicating whether a sequence is ciruclar.
#' @param genome A character vector of length one indicating the genome build.
#'
#' @return a GRanges object
#' @importFrom GenomeInfoDb Seqinfo
#' @export
genome_info <- function(seqnames = NULL, seqlengths = NULL, is_circular = NULL,
                        genome = NULL) {

  seqlengths <- ifelse(length(seqlengths) == 0, NA, seqlengths)
  is_circular <- ifelse(length(is_circular) == 0, NA, is_circular)
  genome <- ifelse(length(genome) == 0, NA, genome)
  seqinfo_res <- Seqinfo(seqnames, seqlengths, is_circular, genome)
  get_genome_info(seqinfo_res)

}

#' Add annotation information as a metadata column to a Ranges
#'
#' @param .data a Ranges object
#' @param ... annotation columns to add to .data, valid names are seqlengths,
#' is_circular, genome, seqnames
#' @export
set_genome_info <- function(.data, ...) { UseMethod("set_genome_info") }

# two options here either use this as a mutate like function so a user
# can construct their own annotations in which case we have name-value pairs
# or we allow to join up the seqinfo part of the ranges object.
# my preference is for the latter, if you want to add in your own columns then
# you should just use mutate
set_genome_info.GenomicRanges <- function(.data, ...) {

  valid_cols <- c("seqnames", "seqlengths", "is_circular", "genome")
  dots <- quos(...)
  dot_names <- unlist(lapply(dots, quo_name))

  if (any(!(dot_names %in% valid_cols))) {
    stop(paste("Invalid selection, can only add",
               paste0(valid_cols,collapse = ","), "."))
  }

  info <- seqinfo(.data)
  seqinfo_df <- data.frame(seqnames = seqnames(info),
                           .seqlengths = seqlengths(info),
                           .is_circular = isCircular(info),
                           .genome = genome(info),
                           stringsAsFactors = FALSE)

  dot_names <- paste0(".", dot_names)
  seqinfo_df <- seqinfo_df[, c("seqnames", dot_names)]
  key <- data.frame(seqnames = seqnames(.data), stringsAsFactors = FALSE)

  annotation_df <- merge(key, seqinfo_df,
                         all.x = TRUE,
                         by = "seqnames", sort = FALSE )
  mcols(.data) <- cbind(mcols(.data), annotation_df[, dot_names, drop = FALSE])

  return(.data)

}

#' Extract annotation as a Ranges object
#' @param .data a Seqinfo object or a GRanges with a seqinfo slot
#' @return A GRanges object
get_genome_info <- function(.data)  UseMethod("get_genome_info")

get_genome_info.GenomicRanges <- function(.data) {
  info <- seqinfo(.data)
  get_genome_info(info)
}


# Similar to as(seqinfo, "GRanges") but adds the information as
# metadata columns.
get_genome_info.Seqinfo <- function(.data) {

  if (any(is.na(seqlengths(.data)))) {
    stop("seqlengths must be non-missing to convert to GRanges")
  }
  gr <- as(.data, "GRanges")
  mcols(gr)$seqlengths <- seqlengths(.data)
  mcols(gr)$is_circular <- isCircular(.data)
  mcols(gr)$genome <- genome(.data)
  gr
}

