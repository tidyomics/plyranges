context("reading/writing bed files")
# tests adapted from rtracklayer

createCorrectGR <- function(seqinfo) {
  ir <- IRanges(c(127471197, 127472364, 127473531, 127474698, 127475865),
                width = 1167)
  space <- factor(rep(c("chr7", "chr9"), c(3, 2)), seqlevels(seqinfo))
  blocks <- S4Vectors::split(IRanges(c(1, 501, 1068, 1, 668, 1, 1, 1),
                                     c(300, 700, 1167, 250, 1167, 1167, 1167, 1167)),
                             rep(seq_len(5), c(3, 2, 1, 1, 1)))
  names(blocks) <- NULL
  correct_gr <- GRanges(space, ir,
                        strand = strand(c("+", "+", "-", "+", "-")),
                        name = c("Pos1", "Pos2", "Neg1", "Pos3", "Neg2"),
                        score = c(0, 2, 0, 5, 5),
                        itemRgb = c("#FF0000", "#FF0000", "#FF0000",
                                    "#FF0000", "#0000FF"),
                        thick = ir, blocks)
  seqinfo(correct_gr) <- seqinfo
  correct_gr
}



test_that("read_bed returns correct GRanges",{

  test_path <- system.file("tests", package = "rtracklayer")

  test_bed <- file.path(test_path, "test.bed")
  correct_gr <- createCorrectGR(Seqinfo(c("chr7", "chr9")))

  test_gr <- read_bed(test_bed)
  # should always return GRanges
  expect_s4_class(test_gr, "GRanges")
  expect_identical(correct_gr, test_gr)

  test_bed_file <- rtracklayer::BEDFile(test_bed)
  test_gr <- read_bed(test_bed_file)
  expect_identical(correct_gr, test_gr)

  test_bed_con <- file(test_bed)
  test_gr <- read_bed(test_bed_con)
  expect_identical(correct_gr, test_gr)

  # check colnames
  subcols <- c("name", "thick")
  # drops strand info
  test_gr <- read_bed(test_bed, col_names = subcols)

  expect_identical(correct_gr %>%
                     mutate(strand = "*") %>%
                     select(subcols),
                   test_gr)
  test_gr <- read_bed(test_bed, col_names = c(subcols, "strand"))
  expect_identical(correct_gr %>%
                     select(subcols),
                   test_gr)

  # check overlaps
  which <- correct_gr[1:2]
  correct_which <- filter_by_overlaps(correct_gr, which)
  test_gr <- read_bed(test_bed, overlap_ranges = which)
  expect_identical(correct_which, test_gr)

  # check that import by genome name works
  if (!require(BSgenome.Hsapiens.UCSC.hg19)) {
    stop("'BSgenome.Hsapiens.UCSC.hg19' must be installed to run tests")
  }

  hg19_seqinfo <- SeqinfoForBSGenome(genome = "hg19")
  correct_gr <- createCorrectGR(hg19_seqinfo)
  test_gr <- read_bed(test_bed, genome_info = "hg19")
  expect_identical(correct_gr, test_gr)

  # if genome_info is a ranges object
  hg19_gr <- get_genome_info(hg19_seqinfo)
  test_gr <- read_bed(test_bed, genome_info = hg19_gr)
  expect_identical(correct_gr, test_gr)

})

test_that("write bed returns correct bed files", {

  correct_gr <- createCorrectGR(Seqinfo(c("chr7", "chr9")))
  ## the 'gsub' is to handle Windows paths (for later coercion to URL)
  dir <- tempdir()
  test_bed_out <- gsub("\\\\", "/", file.path(dir, "test.bed"))
  on.exit(unlink(test_bed_out))

  write_bed(correct_gr, test_bed_out)
  test_gr <- read_bed(test_bed_out)
  expect_identical(correct_gr, test_gr)

  # url input
  test_bed_url <- paste("file://", test_bed_out, sep = "")
  write_bed(correct_gr, test_bed_url)
  test_gr <- read_bed(test_bed_url)
  expect_identical(correct_gr, test_gr)

  # gzipped output
  test_bed_gz <- paste(test_bed_out, ".gz", sep = "")
  write_bed(correct_gr, test_bed_gz)
  test_gr <- read_bed(test_bed_gz)
  expect_identical(correct_gr, test_gr)

  test_bed_gz_url <- paste("file://", test_bed_gz, sep = "")
  write_bed(correct_gr, test_bed_gz_url)
  test_gr <- read_bed(test_bed_gz_url)
  expect_identical(correct_gr, test_gr)

  # tabixed output
  test_bed_bgz <- paste(test_bed_out, ".bgz", sep = "")
  on.exit(unlink(paste(test_bed_bgz, ".tbi", sep = "")))
  write_bed(correct_gr, test_bed_out, index = TRUE)
  which <- correct_gr[1:2]
  correct_which <- filter_by_overlaps(correct_gr, which)
  test_gr <- read_bed(test_bed_bgz, overlap_ranges =  which)
  expect_identical(correct_which, test_gr)

})

test_that("read_narrowpeaks returns correct GRanges", {
  file <- system.file("extdata", "demo.narrowPeak.gz",  package="rtracklayer")
  gr <- read_narrowpeaks(file, genome_info = "hg19")
  expect_equal(length(gr), 6)
  expect_equal(colnames(mcols(gr)),
              c("name","score","signalValue","pValue","qValue","peak"))

  # matches write output
  test_np_out <- file.path(tempdir(), "test.narrowPeak.gz")
  write_narrowpeaks(gr, test_np_out)
  test_gr <- read_narrowpeaks(test_np_out, genome_info = "hg19")
  expect_identical(gr, test_gr)

})
