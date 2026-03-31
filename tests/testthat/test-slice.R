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

test_that("new slice functionality produces expected results", {
  expect_identical(dplyr::slice_head(rng, n=2), rng[1:2])
  expect_identical(dplyr::slice_tail(rng, n=2), rng[9:10])
  expect_identical(dplyr::slice_max(rng, gc), rng[which.max(rng$gc)])
  expect_identical(dplyr::slice_min(rng, gc), rng[which.min(rng$gc)])

  expect_true(length(dplyr::slice_sample(rng, prop=.3)) == 3)
  
  # check that proportion weights look right:
  rng_wts <- rng
  rng_wts$weight <- c(1,1,1,1,1,2,2,2,2,2)
  rep_wts <- replicate(100,
    dplyr::slice_sample(rng_wts, prop=.5, weight_by=weight)$weight
  )
  tab <- table(as.vector(rep_wts))/500
  expect_true(tab[1] < tab[2])

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

  expect_identical(dplyr::slice_head(by_strand) |> ungroup(), rng[c(1,2,5)])
  expect_identical(dplyr::slice_tail(by_strand) |> ungroup(), rng[c(7,9,10)])

  spl <- split(rng$gc, strand(rng))
  target <- dplyr::slice_max(by_strand, gc) %>% ungroup()
  idx <- as.character(strand(target))
  expect_identical(target$gc, sapply(spl[idx], max) %>% unname())
  target <- dplyr::slice_min(by_strand, gc) %>% ungroup()
  idx <- as.character(strand(target))
  expect_identical(target$gc, sapply(spl[idx], min) %>% unname())

  set.seed(123)
  target <- dplyr::slice_sample(by_strand, prop=.5)
  expect_equal(as.numeric(table(strand(target))), c(1,2,1))
  
})