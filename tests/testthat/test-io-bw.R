context("BigWig files")

# tests adapated from rtracklayer
createCorrectGR <- function() {
  ir <- as(IRanges::PartitioningByWidth(rep(300, 9)), "IRanges")
  space <- factor(c(rep("chr2", 5), rep("chr19", 4)), c("chr2", "chr19"))
  score <- seq(-1, 1, length = 9)
  correct_fixed <- GRanges(space, ir, score = score)
  si <- rtracklayer::SeqinfoForBSGenome("hg19")
  seqlengths(correct_fixed) <- seqlengths(si)[levels(space)]
  correct_fixed
}

test_that("reading/ writing bigwig files returns correct GRanges", {
  skip_on_os(os  = "windows")
  test_path <- system.file("tests", package = "rtracklayer")
  test_bw <- file.path(test_path, "test.bw")

  correct_gr <- createCorrectGR()
  test_gr <- read_bigwig(test_bw)

  expect_identical(correct_gr, test_gr)

  # basic writing
  test_bw_out <- file.path(tempdir(), "test_out.bw")
  write_bigwig(correct_gr, test_bw_out)
  on.exit(unlink(test_bw_out))
  test_gr <- read_bigwig(test_bw_out)
  expect_identical(test_gr, correct_gr)

  # bedGraph
  correct_bedgraph <- correct_gr
  width(correct_bedgraph) <- seq(1, 300, length = 9)

  write_bigwig(correct_bedgraph, test_bw_out)
  test_gr <- read_bigwig(test_bw_out)
  expect_identical(test_gr, correct_bedgraph)

  ## overlap ranges
  which <- GRanges(c("chr2", "chr2"), IRanges(c(1, 300), c(400, 1000)))
  correct_which <- filter_by_overlaps(correct_bedgraph, which)
  ranges(correct_which) <- ranges(intersect(correct_which, which))
  test_gr <- read_bigwig(test_bw_out, overlap_ranges = correct_which)
  expect_identical(test_gr, correct_which)

  ## empty overlap
  which <- GRanges()
  correct_which <- filter_by_overlaps(correct_bedgraph, which)
  test_gr <- read_bigwig(test_bw_out, overlap_ranges = which)
  expect_identical(test_gr, correct_which)

})

