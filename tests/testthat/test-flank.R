context("flanks")

test_that("matches HelloRanges", {
  setwd(system.file("unitTests", "data", "flank", package="HelloRanges"))
  gr_a <- read_bed("a.bed")
  left <- flank(gr_a, 5L, ignore.strand = TRUE)
  right <- flank(gr_a, 5L, start = FALSE, ignore.strand = TRUE)

  expect_identical(flank_left(gr_a, 5L), left)
  expect_identical(flank_right(gr_a, 5L), right)

  up <- flank(gr_a, 5L, ignore.strand = FALSE)
  expect_identical(flank_upstream(gr_a, 5L), up)
  down <- flank(gr_a, 5L, start = FALSE, ignore.strand = FALSE)
  expect_identical(flank_downstream(gr_a, 5L), down)

})

test_that("matches GRanges tests", {
  expect_identical(flank_left(GRanges(), 10), GRanges())

  gr <- GRanges(seqnames = c("chr1", "chr2", "chr1", "chrM"),
                ranges = IRanges(21:24, width=10),
                strand = strand(c("+", "-", "*", "-")))

  target <- GRanges(seqnames = seqnames(gr),
                    ranges = IRanges(c(11, 32, 13, 34), width=10),
                    strand = strand(gr))
  current <- flank_upstream(gr, 10L)
  expect_identical(target, current)

  target <- GRanges(seqnames = seqnames(gr),
                    ranges = IRanges(c(31, 12, 33, 14), width=10),
                    strand = strand(gr))
  current <- flank_downstream(gr, 10L)
  expect_identical(target, current)

  target <- GRanges(seqnames = seqnames(gr),
                    ranges = IRanges(c(11, 12, 13, 14), width=10),
                    strand = strand(gr))
  current <- flank_left(gr, 10L)
  expect_identical(target, current)

  target <- GRanges(seqnames = seqnames(gr),
                    ranges = IRanges(c(31, 32, 33, 34), width=10),
                    strand = strand(gr))
  current <- flank_right(gr, 10L)
  expect_identical(target, current)
})
