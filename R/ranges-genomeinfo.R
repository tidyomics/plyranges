# Representing seqinfo information as a Ranges

#' Construct annotation information for Ranges

genome_info <- function(...) {
  valid_names <- c("seqnames", "seqlengths", "is_circular", "genome")
  dots <- rlang::exprs(...)
  if (!any(names(dots) %in% valid_names)) {
    stop("Invalid argument supplied to genome_info.")
  }

  seqinfo_expr <- quo(GenomeInfoDb::Seqinfo(UQS(dots)))
  seqinfo_res <- eval_tidy(seqinfo_expr)
  seqinfo_to_ranges(seqinfo_res)

}

#' Add annotation information as a metadata column to a GRanges
add_genome_info <- function(.data, ...) { UseMethod("add_genome_info") }

#' Extract annotation as a Ranges object
extract_genome_info <- function(.data) { UseMethod("extract_genome_info") }

#' Coerce a Seqinfo object to a GRanges
seqinfo_to_ranges <- function(.data) { UseMethod("seqinfo_to_ranges") }

seqinfo_to_ranges.Seqinfo <- function(.data) {
  if (any(is.na(seqlengths(.data)))) {
    stop("seqlengths must be non-missing to convert to GRanges")
  }
  gr <- GRanges(seqnames = seqnames(.data),
                ranges = IRanges(start = rep(1L, length(.data)),
                                 end = seqlengths(.data)))
  mcols(gr)$seqlengths <- seqlengths(.data)
  mcols(gr)$is_circular <- isCircular(.data)
  mcols(gr)$genome <- genome(.data)
  seqinfo(gr) <- .data
  gr
}

seqinfo_to_ranges.GRanges <- function(.data) {
  info <- seqinfo(.data)
  seqinfo_to_ranges(info)
}
