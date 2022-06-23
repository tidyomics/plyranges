# test-filter-by-overlap.R

test_that("filters by overlap & non-overlap work", {

  subject_gr <- as_granges(data.frame(seqnames = rep("chr1", 3),
                               start = c(1,1, 100),
                               width = rep(100, 3),
                               strand = c("+", "-", "+")))

  query_gr <- as_granges(data.frame(seqnames = c("chr1", "chr2", "chr1"),
                               start = c(20, 20, 150),
                               width = c(10, 10, 10),
                               strand = c("+", "-", "+")))
  
  filter_directed_expected_gr <- as_granges(data.frame(seqnames = rep("chr1", 2),
                               start = c(1, 100),
                               width = rep(100, 2),
                               strand = c("+", "+")))
  
  filter_non_overlap_expected_gr <- GRanges(seqinfo = Seqinfo(seqnames = "chr1"))
  
  filter_non_overlap_directed_expected_gr <- as_granges(
    data.frame(seqnames = "chr1",
                               start = c(1),
                               width = rep(100),
                               strand = c("-")
               )
    )
  
  expect_identical(subject_gr, filter_by_overlaps(subject_gr, query_gr))
  expect_identical(filter_directed_expected_gr, filter_by_overlaps_directed(subject_gr, query_gr))
  
  expect_identical(filter_non_overlap_expected_gr, filter_by_non_overlaps(subject_gr, query_gr))
  expect_identical(filter_non_overlap_directed_expected_gr, filter_by_non_overlaps_directed(subject_gr, query_gr))

})