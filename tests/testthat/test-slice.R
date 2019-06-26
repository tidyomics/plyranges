
set.seed(100)
df <- data.frame(start = 1:10,
                 width = 5,
                 seqnames = "seq1",
                 strand = sample(c("+", "-", "*"), 10, replace = TRUE),
                 gc = runif(10))

rng <- as_granges(df)
by_strand <- group_by(rng, strand)

test_that("standard slicing works", {
  expect_s4_class(dplyr::slice(rng, 1), "GRanges")
  expect_identical(dplyr::slice(rng, 1:n()), rng)
  expect_identical(dplyr::slice(rng, n()), rng[length(rng)])
  expect_identical(dplyr::slice(rng, which.max(gc)), 
                   rng[which.max(rng$gc)])
  expect_identical(dplyr::slice(rng, which.min(gc)), 
                   rng[which.min(rng$gc)])
})

test_that("standard grouping works", {
  expect_s4_class(dplyr::slice(by_strand, 1), "GroupedGenomicRanges")
  expect_equivalent(dplyr::slice(by_strand, n()), 
                   group_by(rng[c(1,6,3)], strand))
  grl_list <- as(split(rng, strand(rng)), "GRangesList")
  index <- as(lapply(grl_list, 
                     function(i) c(which.min(i$gc), which.max(i$gc))),
              "IntegerList")
  correct_gr <- group_by(unname(BiocGenerics::unlist(grl_list[index])),
                         strand)
  expect_equivalent(dplyr::slice(by_strand, which.min(gc), which.max(gc)),
                   correct_gr)
  
})


