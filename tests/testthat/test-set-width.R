context("resizing ranges")




test_that("resizing with anchoring", {
  ir1 <- IRanges(c(2,5,1), c(3,7,3))

  expect_identical(set_width(ir1, 10),
                 IRanges(c(2, 5, 1), width=10))
  expect_identical(set_width(anchor_end(ir1), width = 10),
                 IRanges(c(-6, -2, -6), width=10))
  expect_identical(set_width(anchor_center(ir1), 10),
                 IRanges(c(-2, 1, -3), width=10))


  gr <- GRanges(seqnames = Rle(factor(c("chr1", "chr2", "chr1", "chr3")),
                               c(1, 3, 2, 4)),
                ranges = IRanges(1:10, width = 10:1),
                strand = Rle(strand(c("-", "+", "*", "+", "-")),
                             c(1, 2, 2, 3, 2)))
  # ignores strand by defualt
  gr_w <- set_width(gr, 10)
  expect_identical(rep(10L, length(gr)), width(gr_w))
  expect_identical(c(1L:10L), start(gr_w))
  ir_c <- ranges(set_width(anchor_center(gr), 10))
  expect_identical(ir_c, IRanges(rep(1:5, each=2), width = 10))
  # anchoring at 5p will fix start on + and fix end on negative
  ir_5p <- ranges(set_width(anchor_5p(gr), 10))
  expect_identical(ir_5p,
                 IRanges(c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 1L),
                         width = 10L))
  # anchoring at 3p will fix start on - and fix end on +
  ir_3p <- ranges(set_width(anchor_3p(gr), 10))
  expect_identical(ir_3p,
                   IRanges(c(rep(1L, 8), 9, 10), width = 10L))

})
