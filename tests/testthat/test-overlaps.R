# test-overlaps.R
# adapted from IRanges/GRanges tests
context("overlaps")

make_subject <- function(class = c("GRanges", "IRanges")) {
  class <- match.arg(class)

  switch(class,
         GRanges = GRanges(seqnames = Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
                          ranges = IRanges(1:10, width = 10:1),
                          strand = Rle(strand(c("-", "+", "+", "-", "-", "-")), c(1, 2, 1, 1, 3, 2)),
                          seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep="")),
                          score = 1:10,
                          GC = seq(1, 0, length=10)),
         IRanges = IRanges(c(2, 2, 10), c(2, 3, 12))
  )
}


make_query <- function(class = c("GRanges", "IRanges")) {
  class <- match.arg(class)

  switch(class,
         GRanges = GRanges(seqnames = c("chr1", "chr3", "chr1"),
                           ranges = IRanges(start=c(5,2,1), end=c(10, 7, 5)),
                           strand = c("+", "-", "-"),
                           match = c("no", "one", "two")),
         IRanges = IRanges(c(1, 4, 9), c(5, 7, 10))
  )
}


test_that("group_by_overlaps matches queryHits", {

  subject_gr <- make_subject("GRanges")
  query_gr <- make_query("GRanges")

  # strand is ignored unless method is directed
  hits <- findOverlaps(query_gr, subject_gr, ignore.strand = TRUE)

  expect_equal(group_by_overlaps(query_gr, subject_gr)$query,
               queryHits(hits))

  subject_ir <- make_subject("IRanges")
  query_ir <- make_query("IRanges")

  hits <- findOverlaps(query_ir, subject_ir)

  expect_equal(mcols(group_by_overlaps(query_ir, subject_ir))$query,
               queryHits(hits))

})
