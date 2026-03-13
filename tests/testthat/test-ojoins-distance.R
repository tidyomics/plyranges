# test-overlaps-distance.R
context("overlap joins with distance")

test_that("adding distance works as expected",{

  a <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(c(11, 66, 101), c(15, 70, 110)),
    strand = c("+", "-", "+"),
    a_name = paste0("a", 1:3),
    a_score = 1:3
  )

  b <- GRanges(
    seqnames = "chr1",
    strand = c("-", "+", "-", "+"),
    ranges = IRanges(c(21, 31, 41, 51), c(30, 40, 50, 60)),
    b_name = paste0("b", 1:4),
    b_score = 3:6
  )

  join_overlap_left(a, b, maxgap=5)
  join_overlap_left_directed(a, b, maxgap=15)

  join_nearest(a, b, distance=TRUE)
  join_nearest_downstream(a, b, distance=TRUE)

  join_overlap_left(a, b, maxgap=15, distance=TRUE)
  join_overlap_left_directed(a, b, maxgap=15, distance=TRUE)

})
