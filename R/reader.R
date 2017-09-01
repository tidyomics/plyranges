# readr read_* equivalents
# rtracklayer and Rsamtools cleaner API, they always return
# a GRanges object, perhaps which should be changed to overlap?
# read_bam

# read_bw
#' @importFrom rtracklayer import.bw
read_bigwig <- function(file, ...) {
  import.bw(file, ...)
}

#' @importFrom rtracklayer import.wig
read_wig <- function(file,  genome = NULL, seqinfo = NULL, which = NULL) {
  import.wig(file, genome = genome, seqinfo = NULL, which = NULL)
}

#' @importFrom rtracklayer import.bed
read_bed <- function(file, col_names = NULL, genome = NULL, seqinfo = NULL,
                     which = NULL, col_extra = NULL) {
  import.bed(file, genome = genome, colnames = col_names, genome = genome,
             seqinfo = seqinfo, which = which, extraCols = col_extra)
}

#' @importFrom rtracklayer import.bedGraph
read_bed_graph <- function(file, col_names = NULL, genome = NULL, seqinfo = NULL,
                           which = NULL, col_extra = NULL) {
  import.bedGraph(file, genome = genome, colnames = col_names, genome = genome,
             seqinfo = seqinfo, which = which, extraCols = col_extra)
}

#' @importFrom rtracklayer import.gff
read_gff <- function(file, col_names = NULL, genome = NULL, which = NULL, feature_type = NULL) {
  import.gff(file, colnames = col_names, genome = genome, which = which,
             feature_type = feature_type)
}

#' @importFrom rtracklayer import.gff1
read_gff1 <- function(file, col_names = NULL, genome = NULL, which = NULL, feature_type = NULL) {
  import.gff1(file, colnames = col_names, genome = genome, which = which,
              feature_type = feature_type)
}

#' @importFrom rtracklayer import.gff2
read_gff2 <- function(file, col_names = NULL, genome = NULL, which = NULL, feature_type = NULL) {
  import.gff2(file, colnames = col_names, genome = genome, which = which,
              feature_type = feature_type)
}

#' @importFrom rtracklayer import.gff3
read_gff3 <- function(file, col_names = NULL, genome = NULL, which = NULL, feature_type = NULL) {
  import.gff3(file, colnames = col_names, genome = genome, which = which,
              feature_type = feature_type)
}

