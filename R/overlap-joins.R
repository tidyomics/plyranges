# overlap_joins.R
# join verbs for Ranges-like objects. These are wrapper functions to
# GenomicRanges overlaps methods, and provide a consistent API for
# comapring ranges objects - how should full_ojoin work

# overlap any joins - these target any type joins for ranges
setGeneric("left_ojoin", function(x, y, ...) standardGeneric("left_ojoin"))
setGeneric("right_ojoin", function(x, y, ...) standardGeneric("right_ojoin"))
setGeneric("inner_ojoin", function(x, y, ...) standardGeneric("inner_ojoin"))
setGeneric("full_ojoin", function(x, y, ...) standardGeneric("full_ojoin"))

# convienence functions
# overlap start joins - these target type="start" joins for ranges
setGeneric("start_left_ojoin",
           function(x, y, ...) standardGeneric("start_left_ojoin"))
setGeneric("start_right_ojoin",
           function(x, y, ...) standardGeneric("start_right_ojoin"))
setGeneric("start_inner_ojoin",
           function(x, y, ...) standardGeneric("start_inner_ojoin"))
setGeneric("start_full_ojoin",
           function(x, y, ...) standardGeneric("start_full_ojoin"))

# overlap end joins - these target type="end" joins for ranges
setGeneric("end_left_ojoin",
           function(x, y, ...) standardGeneric("end_left_ojoin"))
setGeneric("end_right_ojoin",
           function(x, y, ...) standardGeneric("end_right_ojoin"))
setGeneric("end_inner_ojoin",
           function(x, y, ...) standardGeneric("end_inner_ojoin"))
setGeneric("end_full_ojoin",
           function(x, y, ...) standardGeneric("end_full_ojoin"))

# overlap equal joins - these target type = "equal" joins for ranges
setGeneric("equal_left_ojoin",
           function(x, y, ...) standardGeneric("equal_left_ojoin"))
setGeneric("equal_right_ojoin",
           function(x, y, ...) standardGeneric("equal_right_ojoin"))
setGeneric("equal_inner_ojoin",
           function(x, y, ...) standardGeneric("equal_inner_ojoin"))
setGeneric("equal_full_ojoin",
           function(x, y, ...) standardGeneric("equal_full_ojoin"))

# overlap start joins - these target type="start" joins for ranges
setGeneric("within_left_ojoin",
           function(x, y, ...) standardGeneric("within_left_ojoin"))
setGeneric("within_right_ojoin",
           function(x, y, ...) standardGeneric("within_right_ojoin"))
setGeneric("within_inner_ojoin",
           function(x, y, ...) standardGeneric("within_inner_ojoin"))
setGeneric("within_full_ojoin",
           function(x, y, ...) standardGeneric("within_full_ojoin"))

#' @name overlap-join-methods
#' @rdname overlap-joins-methods.Rd
#' @param x,y Objects representing ranges
#' @param maxgap,minoverlap The maximimum gap between intervals as an integer
#' greater than or equal to zero. The minimum amount of overlap between intervals
#' as an integer greater than zero, accounting for the maximum gap.
#' @param type the type of overlap
#' @param ignore.strand If set to FALSE account for strand information when
#' joining by overlapping intervals.
#'
#' @seealso \link[GenomicRanges]{setops-methods}, \link[IRanges]{findOverlaps-methods}
#' @importFrom IRanges findOverlaps queryHits subjectHits
setMethod("left_ojoin",
          signature = c("GRanges", "GRanges"),
          function(x, y, maxgap = 0L, minoverlap = 1L,  type = "any", ignore.strand = TRUE) {
            overlaps <- findOverlaps(x, y,
                                      type = type,
                                      maxgap = maxgap,
                                      minoverlap = minoverlap,
                                      ignore.strand = ignore.strand)
            left <- x[queryHits(overlaps), ]
            mcols(left) <- cbind(mcols(left),
                                 mcols(y[subjectHits(overlaps), ]))
            left
          }
)

#' @name overlap-join-methods
#' @rdname overlap-joins-methods.Rd
setMethod("start_left_ojoin", c("GRanges", "GRanges"),
          function(x, y, maxgap = 0L, minoverlap = 1L, ignore.strand = TRUE) {
            left_ojoin(x,y, maxgap, minoverlap, type = "start", ignore.strand = TRUE)
          })

#' @name overlap-join-methods
#' @rdname overlap-joins-methods.Rd
setMethod("end_left_ojoin", c("GRanges", "GRanges"),
          function(x, y, maxgap = 0L, minoverlap = 1L, ignore.strand = TRUE) {
            left_ojoin(x,y, maxgap, minoverlap, type = "end", ignore.strand = TRUE)
          })
#' @name overlap-join-methods
#' @rdname overlap-joins-methods.Rd
setMethod("equal_left_ojoin", c("GRanges", "GRanges"),
          function(x, y, maxgap = 0L, minoverlap = 1L, ignore.strand = TRUE) {
            left_ojoin(x,y, maxgap, minoverlap, type = "equal", ignore.strand = TRUE)
          })

#' @name overlap-join-methods
#' @rdname overlap-joins-methods.Rd
setMethod("within_left_ojoin", c("GRanges", "GRanges"),
          function(x, y, maxgap = 0L, minoverlap = 1L, ignore.strand = TRUE) {
            left_ojoin(x,y, maxgap, minoverlap, type = "within", ignore.strand = TRUE)
          })
#' @name overlap-join-methods
#' @rdname overlap-joins-methods.Rd
#' @importFrom IRanges findOverlapPairs pintersect
setMethod("inner_ojoin",
          signature = c("GRanges", "GRanges"),
          function(x, y, maxgap = 0L, minoverlap = 1L, type = "any", ignore.strand = TRUE) {
            pairs <- findOverlapPairs(x, y,
                                      type = type,
                                      maxgap = maxgap,
                                      minoverlap = minoverlap,
                                      ignore.strand = ignore.strand)
            pintersect(pairs, ignore.strand = ignore.strand)
          }
)

#' @name overlap-join-methods
#' @rdname overlap-joins-methods.Rd
setMethod("start_inner_ojoin", c("GRanges", "GRanges"),
          function(x, y, maxgap = 0L, minoverlap = 1L, ignore.strand = TRUE) {
            inner_ojoin(x,y, maxgap, minoverlap, type = "start", ignore.strand = TRUE)
          })

#' @name overlap-join-methods
#' @rdname overlap-joins-methods.Rd
setMethod("end_inner_ojoin", c("GRanges", "GRanges"),
          function(x, y, maxgap = 0L, minoverlap = 1L, ignore.strand = TRUE) {
            inner_ojoin(x,y, maxgap, minoverlap, type = "end", ignore.strand = TRUE)
          })
#' @name overlap-join-methods
#' @rdname overlap-joins-methods.Rd
setMethod("equal_inner_ojoin", c("GRanges", "GRanges"),
          function(x, y, maxgap = 0L, minoverlap = 1L, ignore.strand = TRUE) {
            inner_ojoin(x,y, maxgap, minoverlap, type = "equal", ignore.strand = TRUE)
          })

#' @name overlap-join-methods
#' @rdname overlap-joins-methods.Rd
setMethod("within_inner_ojoin", c("GRanges", "GRanges"),
          function(x, y, maxgap = 0L, minoverlap = 1L, ignore.strand = TRUE) {
            inner_ojoin(x,y, maxgap, minoverlap, type = "within", ignore.strand = TRUE)
          })

#' @name overlap-join-methods
#' @rdname overlap-joins-methods.Rd
#' @importFrom IRanges findOverlapPairs psetdiff
setMethod("right_ojoin",
          signature = c("GRanges", "GRanges"),
          function(x, y, maxgap = 0L, minoverlap = 1L, type = "any", ignore.strand = TRUE) {
            pairs <- findOverlapPairs(y, x,
                                      type = type,
                                      maxgap = maxgap,
                                      minoverlap = minoverlap,
                                      ignore.strand = ignore.strand)

            right <- y[queryHits(overlaps), ]
            mcols(right) <- cbind(mcols(right),
                                 mcols(x[subjectHits(overlaps), ]))

            right
          }
)

#' @name overlap-join-methods
#' @rdname overlap-joins-methods.Rd
setMethod("start_right_ojoin", c("GRanges", "GRanges"),
          function(x, y, maxgap = 0L, minoverlap = 1L, ignore.strand = TRUE) {
            rigth_ojoin(x,y, maxgap, minoverlap, type = "start", ignore.strand = TRUE)
          })

#' @name overlap-join-methods
#' @rdname overlap-joins-methods.Rd
setMethod("end_right_ojoin", c("GRanges", "GRanges"),
          function(x, y, maxgap = 0L, minoverlap = 1L, ignore.strand = TRUE) {
            right_ojoin(x,y, maxgap, minoverlap, type = "end", ignore.strand = TRUE)
          })
#' @name overlap-join-methods
#' @rdname overlap-joins-methods.Rd
setMethod("equal_right_ojoin", c("GRanges", "GRanges"),
          function(x, y, maxgap = 0L, minoverlap = 1L, ignore.strand = TRUE) {
            right_ojoin(x,y, maxgap, minoverlap, type = "equal", ignore.strand = TRUE)
          })

#' @name overlap-join-methods
#' @rdname overlap-joins-methods.Rd
setMethod("within_right_ojoin", c("GRanges", "GRanges"),
          function(x, y, maxgap = 0L, minoverlap = 1L, ignore.strand = TRUE) {
            right_ojoin(x,y, maxgap, minoverlap, type = "within", ignore.strand = TRUE)
          })
