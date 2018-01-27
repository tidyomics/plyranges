# test-overlaps.R
context("overlap joins")

test_that("intersect join returns correct results",{
  a <- GRanges(seqnames = "chr1",
               ranges = IRanges(c(11,101), c(20, 200)),
               name = c("a1", "a2"),
               strand = c("+", "-"),
               score = c(1,2))
  b <- GRanges(seqnames = "chr1",
               strand = c("+", "-", "+", "-"),
               ranges = IRanges(c(21,91,101,201), c(30,101,110,210)),
               name = paste0("b", 1:4),
               score = 1:4)
  exp <- GRanges("chr1", IRanges(c(11, 101), c(20, 200)), name=c("a1", "a2"),
                 score=c(1, 2), strand=c("+", "-"))

})
