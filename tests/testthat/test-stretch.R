context("stretch ranges")

test_that("matches IRanges", {
  ir1 <- IRanges(c(2,5,1), c(3,7,3))
  ir2 <- ir1
  start(ir2) <- start(ir2) - 10
  end(ir2) <- end(ir2) + 10
  expect_identical(stretch(ir1, 10), ir2)
  ir3 <- ir1
  width(ir3) <- width(ir3) + 10
  expect_identical(stretch(anchor_start(ir1), 10), ir3)
  # will cause negative width
  expect_error(stretch(anchor_end(ir1), 10))
  # anchoring twice will produce same result as one stretch call
  expect_identical(stretch(anchor_end(ir3), -10), ir2)
  # centering
  expect_identical(stretch(anchor_center(ir1), 10),
                   IRanges(c(-8, -4, -8), c(13, 16,12)))
})

test_that("matches GenomicRanges", {
  gr1 <- GRanges(seqnames = Rle(factor(c("chr1", "chr2", "chr1", "chr3")),
                               c(1, 3, 2, 4)),
                ranges = IRanges(1:10, width = 10:1),
                strand = Rle(strand(c("-", "+", "*", "+", "-")),
                             c(1, 2, 2, 3, 2)))
  gr2 <- stretch(gr1, 10)
  expect_equal(width(gr2), width(gr1) + 20)

  gr3 <- gr1
  start(gr3[strand(gr3) == "+"]) <- start(gr3[strand(gr3) == "+"]) - 10
  end(gr3[strand(gr3) == "+"]) <- end(gr3[strand(gr3) == "+"]) + 10
  expect_identical(stretch(anchor_3p(gr1), 10), gr3)

  gr4 <- gr1
  start(gr4[strand(gr4) == "-"]) <- start(gr4[strand(gr4) == "-"]) - 10
  end(gr4[strand(gr4) == "-"]) <- end(gr4[strand(gr4) == "-"]) + 10
  expect_identical(stretch(anchor_5p(gr1), 10), gr4)
})
