context("reading BAM files")

test_that("overlaps works as expected", {
  fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
  olap <- GRanges(seqnames = "seq1", ranges = IRanges(1, width = 100))
  result <- read_bam(fl) %>% filter_by_overlaps(olap)
  expect_s4_class(result, "GRangesDeferred")
  # now has cache filled
  expect_true(has_cache_contents(result))
  expect_identical(c(seq1=1575L, seq2=1584L),
                   seqlengths(get_cache(result)))

})

test_that("no index", {
  fl <- system.file("unitTests", "cases", "ex1_noindex.bam",
                    package="Rsamtools")
  result <- read_bam(fl, index = NULL)
  expect_s4_class(result, "GRangesDeferred")
  expect_false(has_cache_contents(result))
  exp <- DataFrame(n = 3271L)
  expect_identical(exp,  result %>% summarise(n = n()))
})

test_that("selection works with pipes", {
  fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
  exp <- DataFrame(nm_sum = 924L)
  result <- read_bam(fl) %>% select(NM) %>% summarise(nm_sum = sum(NM))
  expect_identical(result, exp)
  result <- read_bam(fl) %>% select(FO) %>% summarise(fo_na = all(is.na(FO)))
  expect_true(result$fo_na)
})

