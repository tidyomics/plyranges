context("stretch ranges")

test_that("matches IRanges", {
  ir1 <- IRanges(c(2,5,1), c(3,7,3))
  ir2 <- ir1
  start(ir2) <- mid(ir2) - 5
  width(ir2) <- 11L
  expect_identical(stretch(ir1, 10), ir2)
  ir3 <- ir1
  width(ir3) <- width(ir3) + 10
  expect_identical(stretch(anchor_start(ir1), 10), ir3)
  # will cause negative width
  expect_error(stretch(anchor_end(ir1), -10))

})

test_that("matches GenomicRanges", {
  # default is to anchor by center
  gr1 <- GRanges(seqnames = Rle(factor(c("chr1", "chr2", "chr1", "chr3")),
                               c(1, 3, 2, 4)),
                ranges = IRanges(1:10, width = 10:1),
                strand = Rle(strand(c("-", "+", "*", "+", "-")),
                             c(1, 2, 2, 3, 2)))
  gr2 <- stretch(gr1, 10)
  expect_equal(mid(gr1), mid(gr2))

  # anchoring by 3p fixes start of - and end of +
  gr3 <- gr1
  pos <- strand(gr3) == "+" | strand(gr3) == "*"
  start(gr3[pos]) <- start(gr3[pos]) - 10
  end(gr3[strand(gr3) == "-"]) <- end(gr3[strand(gr3) == "-"]) + 10
  expect_identical(stretch(anchor_3p(gr1), 10), gr3)

  # anchoring by 3p fixes start of + and end of -
  gr4 <- gr1
  pos <- strand(gr4) == "+" | strand(gr4) == "*"
  start(gr4[strand(gr4) == "-"]) <- start(gr4[strand(gr4) == "-"]) - 10
  end(gr4[pos]) <- end(gr4[pos]) + 10
  expect_identical(stretch(anchor_5p(gr1), 10), gr4)
})
