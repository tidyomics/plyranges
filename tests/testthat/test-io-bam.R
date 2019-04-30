context("reading BAM files")

test_that("overlaps works as expected", {
  fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
  olap <- GRanges(seqnames = "seq1", ranges = IRanges(1, width = 100))
  result <- read_bam(fl) %>% filter_by_overlaps(olap)
  expect_s4_class(result, "DeferredGenomicRanges")
  # now has cache filled
  expect_true(!is_empty_delegate(result))
  # cache has seqlengths methods
  expect_identical(c(seq1=1575L, seq2=1584L), seqlengths(result))
  # can we override cache by selecting a column not present
  result2 <- result %>% select(mapq)
  expect_true("mapq" %in% names(mcols(result2)))
  # but don't override cache if taking out a column
  result3 <- result %>% select(-cigar)
  expect_false("cigar" %in% names(mcols(result3)))
  # filter checks, removing minus strand should result in 1647 Ranges
  result4 <- read_bam(fl) %>% filter(!is_minus_strand)
  expect_true(length(result4) == 1647L)
  # override cache
  result5 <- read_bam(fl) %>% select(seq) %>% filter(!is_unmapped_query)
  expect_true(length(result5) == 3271L)

  # mapq filters
  result6 <- read_bam(fl) %>% filter(mapq > 30)
  expect_true(length(result6) == 3210L)
  
  # proper pairs 
  result7 <- read_bam(fl, paired = TRUE) %>%
    filter(is_proper_pair, !is_duplicate, !is_secondary_alignment, mapq > 30)
  expect_true(length(result7) == 3102L)
  
})

test_that("no index", {
  fl <- system.file("unitTests", "cases", "ex1_noindex.bam",
                    package="Rsamtools")
  result <- read_bam(fl, index = NULL)
  expect_s4_class(result, "DeferredGenomicRanges")
  expect_true(is_empty_delegate(result))
  exp <- DataFrame(n = 3271L)
  expect_identical(exp,  result %>% summarise(n = n()))
})

test_that("selection works with pipes", {
  fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
  exp <- DataFrame(nm_sum = 924L)
  result <- read_bam(fl) %>% 
    select(NM) %>% 
    summarise(nm_sum = sum(NM))
  expect_identical(result, exp)
  
  result <- read_bam(fl) %>% 
    select(FO) %>% 
    summarise(fo_na = all(is.na(FO)))
  expect_true(result$fo_na)
})