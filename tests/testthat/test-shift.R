# test-shift.R
# Adapted from HelloRanges tests.

context("shift")
stopifnot(requireNamespace("HelloRanges"))
oldwd <- getwd()
setwd(system.file("unitTests", "data", "shift", package="HelloRanges"))

genome <- read.delim("tiny.genome", header = FALSE)
genome <- genome_info(genome = "test",
                      seqnames = as.character(genome$V1),
                      seqlengths = genome$V2)

test_that("compatible with GenomicRanges shift:", {
  gr <- read_bed("a.bed", genome_info = genome)
  exp <- shift(gr, 5L)
  expect_identical(exp, shift_right(gr, 5L))
  exp <- gr
  exp[strand(exp) == "-"] <- shift(exp[strand(exp) == "-"], 5L)
  expect_identical(exp, shift_upstream(anchor_3p(gr), 5L))
  exp <- gr
  exp[strand(exp) == "+"] <- shift(exp[strand(exp) == "+"], 5L)
  expect_identical(exp, shift_downstream(anchor_5p(gr), 5L))
})

setwd(oldwd)
