set.seed(1999)
df <- data.frame(start = 1:10,
                 width = 5,
                 seqnames = "seq1",
                 strand = sample(c("+", "-", "*"), 10, replace = TRUE),
                 gc = runif(10))
rng <- as_granges(df)

test_that("rangewise slice produces expected results", {
  expect_identical(dplyr::slice(rng, 1:2), rng[1:2])
  expect_identical(dplyr::slice(rng, -n()), rng[1:9])
  expect_identical(dplyr::slice(rng, -5:-n()), rng[1:4])
  expect_error(dplyr::slice(rng, n()+1L))
  expect_error(dplyr::slice(rng, gc > 0.5))
})



by_strand <- group_by(rng, strand)

test_that("groupwise slice produces expected results", {

  target <- dplyr::slice(by_strand, n()) %>% ungroup()
  exp <- rng[c(7,9,10)]
  expect_identical(target, exp)
  
  
  target <- dplyr::slice(by_strand, which.max(gc))
  exp <- filter(by_strand, gc == max(gc))
  expect_identical(target, exp)
  
  target <- dplyr::slice(by_strand, 1:4)
  expect_identical(target, by_strand)
  
})