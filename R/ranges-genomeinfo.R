# Representing seqinfo information as a Ranges

#' @name ranges-info
#' @rdname ranges-info
#' @export
genome_info <- function(genome = NULL, seqnames = NULL, seqlengths = NULL, is_circular = NULL) {

  if (length(seqlengths) == 0L) seqlengths <- NA
  if (length(is_circular) == 0L) is_circular <- NA
  if (length(genome) == 0L) genome <- NA

  seqinfo_res <- Seqinfo(seqnames, seqlengths, is_circular, genome)
  get_genome_info(seqinfo_res)

}

#' Construct annotation information
#'
#' @param .data A Ranges object to annotate or retrieve an annotation for.
#' @param genome A character vector of length one indicating the genome build.
#' If this is the only argument supplied, the build information will be
#' retrieved from UCSC database.
#' @param seqnames A character vector containing the name of sequences.
#' @param seqlengths An optional integer vector containg the lengths of sequences.
#' @param is_circular An optional logical vector indicating whether a sequence is ciruclar.
#'
#' @return a GRanges object containing annotations. To retrieve the annotations
#' as a Ranges object use \code{get_genome_info}.
#'
#' @description To construct annotations by supplying annotation information
#' use \code{genome_info}. This function allows you to get UCSC build information
#' via \link[GenomeInfoDb]{fetchExtendedChromInfoFromUCSC}. To add
#' annotations to an existing Ranges object use \code{set_genome_info}. To retrieve
#' an annotation as a Ranges object use \code{get_genome_info}.
#'
#' @importFrom GenomeInfoDb Seqinfo seqnames seqlengths isCircular genome seqinfo
#' @seealso \link[GenomeInfoDb]{Seqinfo-class} \link[GenomeInfoDb]{fetchExtendedChromInfoFromUCSC}
#' @examples
#' if (interactive()) {
#'  # requires internet connection
#'  genome_info(genome = "hg38")
#' }
#'
#' x <- genome_info(genome = "toy",
#'                  seqnames = letters[1:4],
#'                  seqlengths = c(100, 300, 15, 600),
#'                  is_circular = c(NA, FALSE, FALSE, TRUE))
#'x
#'
#' rng <- as_granges(data.frame(seqnames = "a", start = 30:50, width = 10))
#' rng
#' rng <- set_genome_info(rng,
#'                        genome = "toy",
#'                        seqnames = letters[1:4],
#'                        seqlengths = c(100, 300, 15, 600),
#'                        is_circular = c(NA, FALSE, FALSE, TRUE))
#' get_genome_info(rng)
#'
#'
#' @rdname ranges-info
#' @export
set_genome_info <- function(.data, genome = NULL, seqnames = NULL,
                            seqlengths = NULL, is_circular = NULL) {



  if (!is.null(genome)) {
    GenomeInfoDb::genome(.data) <- genome
  }

  if (!is.null(seqlengths)) {
    GenomeInfoDb::seqlengths(.data) <- genome
  }

  if (!is.null(seqnames)) {
    if (any(!(seqnames %in% seqnames(.data)))) {
      stop("Provide seqnames do not match seqnames(.data).", call. = FALSE)
    }
    GenomeInfoDb::seqnames(.data) <- seqnames
  }

  if (!is.null(is_circular)) {
    GenomeInfoDb::isCircular(.data) <- is_circular
  }

  return(.data)
}



#' @export
#' @importFrom methods hasMethod
#' @rdname ranges-info
get_genome_info <- function(.data)  UseMethod("get_genome_info")

#' @export
get_genome_info.default <- function(.data) {
  if (hasMethod("seqinfo", class(.data))) {
    return(get_genome_info.Seqinfo(seqinfo(.data)))
  } else {
    stop(paste("Unable to retrieve genome annotation from object of class:",
         class(.data)), call. = FALSE)
  }
}

#' @export
get_genome_info.Seqinfo <- function(.data) {

  if (any(is.na(seqlengths(.data)))) {
    stop("seqlengths must be non-missing to convert to GRanges")
  }
  as(.data, "GRanges")
}

