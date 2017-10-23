# test-coverage.R
context("coverage")
stopifnot(requireNamespace("HelloRanges"))

oldwd <- getwd()
setwd(system.file("unitTests", "data", "genomecov", package="HelloRanges"))

genome <- read.delim("test.genome", header = FALSE)

genome <- genome_info(genome = "test",
                      seqnames = as.character(genome$V1),
                      seqlengths = genome$V2)

test_that("Compatible with GenomicRanges", {
  y <- read_bed("y.bed", genome_info = genome)
  cov_gr <- as(coverage(granges(y)), "GRanges")
  exp <- cov_gr[score(cov_gr) > 0]
  cov_plyr <- y %>% set_coverage() %>% filter(score > 0)
  expect_identical(exp, cov_plyr)

})

setwd(oldwd)
