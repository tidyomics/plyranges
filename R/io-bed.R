# io-bed.R
norm_args_reader <- function(genome_info) {

  if (is(genome_info, "GenomicRanges")) {
    seq_info <- seqinfo(genome_info)
    genome_info <- NA
  } else if (is.character(genome_info)) {
    seq_info <- NULL
  } else {
    genome_info <- NA
    seq_info <- NULL
  }

  return(list(genome_info = genome_info, seq_info = seq_info))
}

#' Read a BED or BEDGraph file
#'
#' @param file A path to a file or a connection.
#' @param col_names An optional character vector for including additional
#' columns in \code{file} that are not part of the BED/narrowPeaks specification.
#' @param genome_info An optional character string or a Ranges object
#' that contains information about the genome build. For example the USSC identifier
#'"hg19" will add build information to the returned GRanges.
#' @param overlap_ranges An optional Ranges object. Only the intervals in the file
#' that overlap the Ranges will be returned.
#'
#' @description This is a lightweight wrapper to the import family
#' of functions defined in \pkg{rtracklayer}.
#'
#' @description Read common interval based formats as GRanges.
#'
#' @details This is a lightweight wrapper to the import family
#' of functions defined in \pkg{rtracklayer}.
#' The \code{read_narrowpeaks} function parses the ENCODE narrowPeak BED format (see
#' \url{https://genome.ucsc.edu/FAQ/FAQformat.html#format12} for details.). As
#' such the parser expects four additional columns called (corresponding to
#' the narrowPeaks spec):
#' \itemize{
#'   \item signalValue
#'   \item pValue
#'   \item qValue
#'   \item peak
#' }
#'
#' @return A GRanges object
#'
#' @importFrom rtracklayer import.bed
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom methods is
#' @seealso \code{\link[rtracklayer]{BEDFile}}
#'
#' @examples
#'
#' test_path <- system.file("tests", package = "rtracklayer")
#  bed_file <- file.path(test_path, "test.bed")
#' gr <- read_bed(bed_file)
#' gr
#' gr <- read_bed(bed_file, genome_info = "hg19")
#' gr
#' olap <-  as_granges(data.frame(seqnames = "chr7", start = 1, end = 127473000))
#' gr <- read_bed(bed_file,
#'               overlap_ranges = olap)
#' # bedGraph
#' bg_file <- file.path(test_path, "test.bedGraph")
#' gr <- read_bedgraph(bg_file)
#' gr
#' # narrowpeaks
#' np_file <- system.file("extdata", "demo.narrowPeak.gz",  package="rtracklayer")
#' gr <- read_narrowpeaks(np_file, genome_info = "hg19")
#' gr
#'
#' @export
#' @rdname io-bed-read
read_bed <- function(file, col_names = NULL, genome_info = NULL,
                     overlap_ranges = NULL) {
  # check genome_info input
  args <- norm_args_reader(genome_info)

  import.bed(file, colnames = col_names,
             genome = args$genome_info,
             seqinfo = args$seq_info,
             which = overlap_ranges,
             trackLine = FALSE)
}


#' Write a BED or BEDGraph file
#'
#' @param x A GRanges object
#' @param file File name, URL or connection specifying a file to write x to.
#'             Compressed files with extensions such as '.gz' are handled
#'             automatically. If you want to index the file with tabix use the
#'             \code{index} argument.
#' @param index Compress and index the output file
#'              with bgzf and tabix (default = FALSE). Note that tabix indexing will sort the
#'              data by chromosome and start.
#'
#' @description This is a lightweight wrapper to the export family
#' of functions defined in \pkg{rtracklayer}.
#'
#' @importFrom rtracklayer export.bed
#' @seealso \link[rtracklayer]{BEDFile} \link[rtracklayer]{BEDGraphFile}
#' @export
#' @rdname io-bed-write
write_bed <- function(x, file, index = FALSE) {
  export.bed(x, file, index = index)
}

#' @rdname io-bed-read
#' @export
#' @importFrom rtracklayer import.bedGraph
read_bed_graph <- function(file, col_names = NULL, genome_info = NULL,
                           overlap_ranges = NULL) {

  args <- norm_args_reader(genome_info)

  import.bedGraph(file, colnames = col_names,
                  genome = args$genome_info,
                  seqinfo = args$seq_info,
                  which = overlap_ranges,
                  trackLine = FALSE)
}


#' @importFrom rtracklayer export.bedGraph
#' @export
#' @rdname io-bed-write
write_bed_graph <- function(x, file, index = FALSE) {
  export.bedGraph(x, file, index = index)
}


#' @export
#' @rdname io-bed-read
read_narrowpeaks <- function(file,
                             col_names = NULL,
                             genome_info = NULL,
                             overlap_ranges = NULL) {

  args <- norm_args_reader(genome_info)

  extra_cols <- c(signalValue="numeric",
                  pValue="numeric",
                  qValue="numeric",
                  peak="integer")

  import.bed(file,
             colnames = col_names,
             extraCols = extra_cols,
             genome = args$genome_info,
             seqinfo = args$seq_info,
             which = overlap_ranges)
}


#' @export
#' @rdname io-bed-write
write_narrowpeaks <- function(x, file) {
  # at the moment just use write.table
  extra_cols <- c("signalValue", "pValue", "qValue", "peak")
  valid_np <- all(extra_cols %in% names(mcols(x)))
  if (!valid_np) {
    stop(paste("For a valid narrowPeaks there must be columns called:",
               paste(extra_cols, collapse = ","), "in x."),
         call. = FALSE)
  }

  # bed files are 0-based
  rng <- DataFrame(chr = seqnames(x),
                   start = start(x) - 1L,
                   end = end(x),
                   strand = gsub("\\*", "\\.", as.character(strand(x))))

  np_df <- mcols(x)[, extra_cols]

  np_name <- try(mcols(x)[, "name", drop = FALSE], silent = TRUE)
  if (is(np_name, "try-error")) {
    np_name <- DataFrame(name = rep(".", nrow(rng)))
  }

  np_score <- try(mcols(x)[, "score", drop = FALSE], silent = TRUE)
  if (is(np_score, "try-error")) {
    np_score <- DataFrame(score = rep(0L, nrow(rng)))
  }
  np_col_order <- c("chr", "start", "end", "name", "score", "strand",
                    extra_cols)

  np_df <- cbind(rng, np_name, np_score, np_df)[, np_col_order]

  utils::write.table(np_df,
                     file,
                     sep = "\t",
                     row.names = FALSE,
                     col.names = FALSE,
                     na = ".",
                     quote = FALSE)
}
