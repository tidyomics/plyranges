# test-coverage.R
# adapted from HelloRanges tests
context("coverage")

test_that("produces same result as IRanges",
          {
            ir <- IRanges(c(1, 8, 14, 15, 19, 34, 40),
                          width = c(12, 6, 6, 15, 6, 2, 7))
            cvg <- compute_coverage(ir)
            expect_identical(mcols(cvg)$score,
                             c(1L, 2L, 1L, 2L, 3L, 2L, 1L, 0L, 1L, 0L, 1L))
            expect_identical(width(cvg),
                             c(7L, 5L, 2L, 4L, 1L, 5L, 5L, 4L, 2L, 4L, 7L))

            ir <- IRanges(start=c(-2L, 6L, 9L, -4L, 1L, 0L, -6L, 10L),
                          width=c( 5L, 0L, 6L,  1L, 4L, 3L,  2L,  3L))
            cvg <- compute_coverage(ir)
            expect_identical(mcols(cvg)$score,
                             c(3L, 1L, 0L, 1L, 2L, 1L))
            expect_identical(width(cvg), c(2L, 2L, 4L, 1L, 3L, 2L))

            cvg <- ir %>% shift_right(7L) %>% compute_coverage()
            expect_identical(mcols(cvg)$score,
                             c(1L, 0L, 1L, 2L, 3L, 1L, 0L, 1L, 2L, 1L))
            expect_identical(width(cvg),
                             c(3L, 1L, 2L, 1L, 2L, 2L, 4L, 1L, 3L, 2L))

          })

test_that("produces same results as GRanges", {
  gr <- GRanges(seqnames = Rle(factor(c("chr1", "chr2", "chr1", "chr3")),
                               c(1, 3, 2, 4)),
                ranges = IRanges(1:10, width = 10:1))
  target <- IRanges::RleList(chr1=Rle(1:3, c(4, 1, 5)),
                    chr2=Rle(0:3, c(1, 1, 1, 7)),
                    chr3=Rle(0:4, c(6, 1, 1, 1, 1)),
                    compress=FALSE)
  cvg <- compute_coverage(gr)

  expect_identical(cvg$score, unlist(runValue(target), use.names = FALSE))
  expect_identical(width(cvg), unlist(runLength(target), use.names = FALSE))

  target <- RleList(chr1=Rle(1:3, c(4, 1, 5)),
                    chr2=Rle(c(0:3, 0L), c(1, 1, 1, 7, 10)),
                    chr3=Rle(c(0:4, 0L), c(6, 1, 1, 1, 1, 20)),
                    compress=FALSE)

  cvg <- compute_coverage(gr, width = c(chr1 = 10, chr2 = 20, chr3 = 30))

  expect_identical(cvg$score, unlist(runValue(target), use.names = FALSE))
  expect_identical(width(cvg), unlist(runLength(target), use.names = FALSE))

})
test_that("Compatible with GenomicRanges", {
  skip_if_not(
    requireNamespace("HelloRanges", quietly = TRUE),
    message = "'HelloRanges' must be installed to run coverage unit tests."
  )
  oldwd <- getwd()
  setwd(system.file("unitTests", "data", "genomecov", package="HelloRanges"))

  genome <- read.delim("test.genome", header = FALSE)

  genome <- genome_info(genome = "test",
                        seqnames = as.character(genome$V1),
                        seqlengths = genome$V2)
  y <- read_bed("y.bed", genome_info = genome)
  cov_gr <- as(coverage(granges(y)), "GRanges")
  exp <- cov_gr[score(cov_gr) > 0]
  cov_plyr <- y %>% compute_coverage() %>% filter(score > 0)
  expect_identical(exp, cov_plyr)
  setwd(oldwd)

})

