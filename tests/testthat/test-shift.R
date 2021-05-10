# test-shift.R
# Adapted from HelloRanges tests.
context("shift")
skip_if_not(
  requireNamespace("HelloRanges", quietly = TRUE),
  message = "'HelloRanges' must be installed to run shift unit tests."
)
oldwd <- getwd()
setwd(system.file("unitTests", "data", "shift", package="HelloRanges"))

genome <- read.delim("tiny.genome", header = FALSE)
genome <- genome_info(genome = "test",
                      seqnames = as.character(genome$V1),
                      seqlengths = genome$V2)

gr <- read_bed("a.bed", genome_info = genome)

test_that("compatible with GenomicRanges shift:", {
  exp <- shift(gr, 5L)
  expect_identical(exp, shift_right(gr, 5L))
  expect_identical(shift(gr, -5L), shift_left(gr, 5L))
  exp <- gr
  exp[strand(exp) == "-"] <- shift(exp[strand(exp) == "-"], 5L)
  exp[strand(exp) == "+"] <- shift(exp[strand(exp) == "+"], -5L)
  expect_identical(exp, shift_upstream(gr, 5L))
  exp <- gr
  exp[strand(exp) == "+"] <- shift(exp[strand(exp) == "+"], 5L)
  exp[strand(exp) == "-"] <- shift(exp[strand(exp) == "-"], -5L)
  expect_identical(exp, shift_downstream(gr, 5L))

})


test_that("shift inputs", {
  # shift must be postive
  expect_error(shift_left(gr, shift = -5))
  # shift must have same length as input ranges
  expect_error(shift_downstream(gr, 1:3))
  # see issue https://github.com/sa-lee/plyranges/issues/73
  # vector input
  exp <- shift(gr, c(5L, 3L))
  expect_identical(shift_right(gr, c(5L, 3L)), exp)
  exp <- shift(gr, c(-5L, -3L))
  expect_identical(shift_left(gr, c(5L, 3L)), exp)
  expect_identical(shift_upstream(gr, c(5,5)), shift_upstream(gr, 5))

  # shift vector is parallel to ranges
  exp <- gr
  exp[strand(exp) == "-"] <- shift(exp[strand(exp) == "-"], 3L)
  exp[strand(exp) == "+"] <- shift(exp[strand(exp) == "+"], -5L)
  expect_identical(exp, shift_upstream(gr, c(5L, 3L)))
  exp <- gr
  exp[strand(exp) == "+"] <- shift(exp[strand(exp) == "+"], 5L)
  exp[strand(exp) == "-"] <- shift(exp[strand(exp) == "-"], -3L)
  expect_identical(exp, shift_downstream(gr, c(5L, 3L)))

  # when there's unstranded ranges
  gr <- GRanges(seqnames="chr1",
                strand=c("+", "-", "+", "*"),
                ranges=IRanges(start=c(3, 8, 11, 5), width=2))
  expect_identical(shift_downstream(gr, c(3, 3, 3, 3)),
                   shift_downstream(gr, 3))

  exp <- gr
  exp[strand(exp) == "-"] <- shift(exp[strand(exp) == "-"], 3)
  exp[strand(exp) %in% c("*", "+")] <- shift(exp[strand(exp) %in% c("*", "+")],
                                             c(-0,-1,-2))
  expect_identical(exp, shift_upstream(gr, c(0,3,1,2)))

  exp <- gr
  exp[strand(exp) == "-"] <- shift(exp[strand(exp) == "-"], -3)
  exp[strand(exp) %in% c("*", "+")] <- shift(exp[strand(exp) %in% c("*", "+")],
                                             c(0,1,2))
  expect_identical(exp, shift_downstream(gr, c(0,3,1,2)))
})

setwd(oldwd)
