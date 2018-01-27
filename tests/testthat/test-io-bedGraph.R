context("reading/writing bedGraph files")
# tests adapted from rtracklayer

createCorrectGR <- function(seqinfo) {
  ir <- shift(IRanges(end = seq(300, 2700, by = 300), width = 300), 59102000)
  score <- seq(-1, 1, by = 0.25)
  space <- factor(c(rep("chr19", 6), "chr17", rep("chr18", 2)),
                  seqlevels(seqinfo))
  correct_gr <- GRanges(space, ir, score = score)
  if (!any(is.na(genome(seqinfo))))
    genome(correct_gr) <- unname(genome(seqinfo)[1])
  seqinfo(correct_gr) <- seqinfo
  correct_gr
}

test_that("read_bed_graph returns correct GRanges", {
  test_path <- system.file("tests", package = "rtracklayer")
  test_bg <- file.path(test_path, "test.bedGraph")
  correct_gr <- createCorrectGR(Seqinfo(c("chr19", "chr17", "chr18")))

  test_gr <- read_bed_graph(test_bg)
  expect_identical(test_gr, correct_gr)

  test_bg_file <- rtracklayer::BEDGraphFile(test_bg)
  test_gr <- read_bed_graph(test_bg)
  expect_identical(test_gr, correct_gr)
  test_bg_con <- file(test_bg)
  test_gr <- read_bed_graph(test_bg_con)
  expect_identical(test_gr, correct_gr)


  # check overlaps
  correct_which <- filter_by_overlaps(correct_gr, correct_gr[3:4])
  test_gr <- read_bed_graph(test_bg, overlap_ranges = correct_gr[3:4])
  expect_identical(correct_which, test_gr)

  # test genome_info
  if (!require(BSgenome.Hsapiens.UCSC.hg19)) {
    stop("'BSgenome.Hsapiens.UCSC.hg19' must be installed to run tests")
  }

  hg19_seqinfo <- SeqinfoForBSGenome("hg19")
  correct_genome <- createCorrectGR(hg19_seqinfo)
  test_gr <- read_bed_graph(test_bg, genome_info = "hg19")
  expect_identical(correct_genome, test_gr)

  hg19_gr <- get_genome_info(hg19_seqinfo)
  test_gr <- read_bed_graph(test_bg, genome_info = hg19_gr)
  expect_identical(correct_genome, test_gr)

})

test_that("writing bedGraph files works", {

  correct_gr <- createCorrectGR(Seqinfo(c("chr19", "chr17", "chr18")))
  test_bg_out <- file.path(tempdir(), "test.bedGraph")
  on.exit(unlink(test_bg_out))

  write_bed_graph(correct_gr, test_bg_out)

  test_gr <- read_bed_graph(test_bg_out)
  expect_identical(correct_gr, test_gr)

  test_foo_out <- file.path(tempdir(), "test.foo")
  write_bed_graph(correct_gr, test_foo_out)
  on.exit(unlink(test_foo_out))
  test_gr <- read_bed_graph(test_bg_out)
  expect_identical(correct_gr, test_gr)
  test_bg_out_file <- rtracklayer::BEDGraphFile(test_bg_out)
  write_bed_graph(correct_gr, test_bg_out_file)
  test_gr <- read_bed_graph(test_bg_out_file)
  expect_identical(correct_gr, test_gr)

  # gzipped output
  test_bg_gz <- paste(test_bg_out, ".gz", sep = "")
  on.exit(unlink(test_bg_gz))
  write_bed_graph(correct_gr, test_bg_gz)
  test_gr <- read_bed_graph(test_bg_gz)
  expect_identical(correct_gr, test_gr)

})
