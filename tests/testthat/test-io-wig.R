context("reading/writing WIG files")

createCorrectGR <- function(seqinfo) {
  score <- c(seq(10, 20, by = 2.5), seq(17.5, 10, by = -2.5),
             seq(1000, 100, by = -100))
  start <- c(59104701 + cumsum(c(0, 200, 500, 200, 300, 180, 220, 390, 1180)),
             59108021 + cumsum(c(0, rep(300, 9))))
  width <- c(rep(1, 9), rep(200, 10))
  space <- factor(c(rep("chr19", 9), rep("chr18", 10)), seqlevels(seqinfo))
  correct_gr <- GRanges(space, IRanges(start, width = width), score = score)
  if (!any(is.na(genome(seqinfo))))
    genome(correct_gr) <- unname(genome(seqinfo)[1])
  seqinfo(correct_gr) <- seqinfo
  correct_gr
}

test_that("reading WIG files", {
  test_path <- system.file("tests", package = "rtracklayer")
  test_wig <- file.path(test_path, "step.wig")

  correct_gr <- createCorrectGR(Seqinfo(c("chr19", "chr18")))
  test_gr <- read_wig(test_wig)
  expect_identical(correct_gr, test_gr)
  test_wig_file <- WIGFile(test_wig)
  test_gr <- read_wig(test_wig_file)
  expect_identical(correct_gr, test_gr)

  test_wig_con <- file(test_wig)
  test_gr <- read_wig(test_wig_con)
  expect_identical(correct_gr, test_gr)


  ## genome_info checks
  hg19_seqinfo <- SeqinfoForBSGenome("hg19")
  correct_genome <- createCorrectGR(hg19_seqinfo)
  test_gr <- read_wig(test_wig, genome_info = "hg19")
  expect_identical(correct_genome, test_gr)

  hg19_gr <- get_genome_info(hg19_seqinfo)
  test_gr <- read_wig(test_wig, genome_info = hg19_gr)

  expect_identical(correct_genome, test_gr)

  # overlaps
  which <- correct_gr[3:4]
  correct_which <- filter_by_overlaps(correct_gr, which)
  test_gr <- read_wig(test_wig, overlap_ranges = which)
  expect_identical(correct_which, test_gr)

})

test_that("writing wig files", {

  correct_gr <- createCorrectGR(Seqinfo(c("chr19", "chr18")))

  test_wig_out <- file.path(tempdir(), "test.wig")
  on.exit(unlink(test_wig_out))
  write_wig(correct_gr, test_wig_out)
  test_gr <-read_wig(test_wig_out)
  expect_identical(correct_gr, test_gr)

  test_foo_out <- file.path(tempdir(), "test.foo")
  write_wig(correct_gr, test_foo_out)
  on.exit(unlink(test_foo_out))
  test_gr <- read_wig(test_wig_out)
  expect_identical(test_gr, correct_gr)

  test_wig_out_file <- rtracklayer::WIGFile(test_wig_out)
  write_wig(correct_gr, test_wig_out_file)
  test_gr <- read_wig(test_wig_out_file)
  expect_identical(test_gr, correct_gr)

  test_wig_gz <- paste(test_wig_out, ".gz", sep = "")
  on.exit(unlink(test_wig_gz))
  write_wig(correct_gr, test_wig_gz)
  test_gr <- read_wig(test_wig_gz)
  expect_identical(test_gr, correct_gr)

})
