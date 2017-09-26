# ranges-joins
# need to include a by argument here
ojoin_inner <- function(x, y, maxgap, minoverlap) UseMethod("ojoin_inner")

ojoin_inner.Ranges <- function(x, y, maxgap = 0L, minoverlap = 1L) {
  pairs <- findOverlapPairs(x, y,
                            type = "any",
                            maxgap = maxgap,
                            minoverlap = minoverlap)
  pintersect(pairs)
}


ojoin_inner.GenomicRanges <- function(x, y, maxgap = 0L, minoverlap = 1L) {
  pairs <- findOverlapPairs(x, y,
                            type = "any",
                            maxgap = maxgap,
                            minoverlap = minoverlap,
                            ignore.strand = TRUE)
  pintersect(pairs, ignore.strand = TRUE)

}


ojoin_left <- function(x, y, maxgap, minoverlap) UseMethod("ojoin_left")

ojoin_left.GenomicRanges <- function(x, y, maxgap = 0L, minoverlap = 1L) {
  overlaps <- findOverlaps(x, y,
                           type = "any",
                           maxgap = maxgap,
                           minoverlap = minoverlap)
  left <- x[queryHits(overlaps), ]
  mcols(left) <- cbind(mcols(left),
                       mcols(y[subjectHits(overlaps), ]))
  left
}

ojoin_left.GenomicRanges <- function(x, y, maxgap = 0L, minoverlap = 1L) {
  overlaps <- findOverlaps(x, y,
                           type = "any",
                           maxgap = maxgap,
                           minoverlap = minoverlap,
                           ignore.strand = TRUE)
  left <- x[queryHits(overlaps), ]
  mcols(left) <- cbind(mcols(left),
                       mcols(y[subjectHits(overlaps), ]))
  left
}
