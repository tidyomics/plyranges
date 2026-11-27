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

test_that("slice_sample, slice_min, and slice_max work with core columns and expressions", {
  gr_var <- as_granges(data.frame(
    seqnames = c("chr1", "chr1", "chr2", "chr2"),
    start = c(10, 1, 20, 5),
    width = c(1, 1000, 1, 1000),
    gc = c(0.1, 0.9, 0.2, 0.8)
  ))

  set.seed(42)
  draws <- replicate(50, width(dplyr::slice_sample(gr_var, n = 1, weight_by = width)))
  expect_true(mean(draws == 1000L) > 0.9)

  expect_identical(dplyr::slice_max(gr_var, order_by = width, n = 1), gr_var[c(2, 4)])
  expect_identical(dplyr::slice_min(gr_var, order_by = start, n = 1), gr_var[2])
  expect_identical(dplyr::slice_max(gr_var, order_by = end, n = 1), gr_var[4])
  expect_identical(dplyr::slice_max(gr_var, order_by = width * gc, n = 1), gr_var[2])

  by_seq <- group_by(gr_var, seqnames)
  expect_identical(dplyr::slice_max(by_seq, order_by = width, n = 1) |> ungroup(), gr_var[c(2, 4)])
  expect_identical(dplyr::slice_min(by_seq, order_by = start, n = 1) |> ungroup(), gr_var[c(2, 4)])

  set.seed(42)
  grp_draws <- replicate(
    50,
    width(dplyr::slice_sample(by_seq, n = 1, weight_by = width) |> ungroup())
  )
  expect_true(mean(grp_draws == 1000L) > 0.9)

  ir_var <- as_iranges(data.frame(start = c(10, 1, 5), width = c(1, 1000, 1)))
  expect_identical(dplyr::slice_max(ir_var, order_by = width, n = 1), ir_var[2])
  expect_identical(dplyr::slice_min(ir_var, order_by = start, n = 1), ir_var[2])
  expect_length(dplyr::slice_sample(ir_var, n = 1, weight_by = width), 1L)
})

test_that("slice_sample, slice_min, and slice_max work with Rle columns", {
  gr_rle <- GRanges(
    seqnames = c("chr1", "chr1", "chr2", "chr2"),
    ranges = IRanges(start = c(10, 1, 20, 5), width = 1),
    score = S4Vectors::Rle(c(1, 5, 2, 8)),
    w = S4Vectors::Rle(c(0, 1, 0, 1))
  )

  expect_identical(dplyr::slice_max(gr_rle, order_by = score, n = 1), gr_rle[4])
  expect_identical(dplyr::slice_min(gr_rle, order_by = score, n = 2), gr_rle[c(1, 3)])
  expect_identical(dplyr::slice_max(gr_rle, order_by = score * width, n = 1), gr_rle[4])

  by_seq <- group_by(gr_rle, seqnames)
  expect_identical(dplyr::slice_max(by_seq, order_by = score, n = 1) |> ungroup(), gr_rle[c(2, 4)])
  expect_identical(dplyr::slice_min(by_seq, order_by = score, n = 1) |> ungroup(), gr_rle[c(1, 3)])

  # ranges with zero weight are never drawn
  expect_identical(sort(dplyr::slice_sample(gr_rle, n = 2, weight_by = w)), sort(gr_rle[c(2, 4)]))
  expect_identical(dplyr::slice_sample(by_seq, n = 1, weight_by = w) |> ungroup(), gr_rle[c(2, 4)])

  # seqnames and strand are factor-Rle columns, not numeric ones
  expect_error(dplyr::slice_max(gr_rle, order_by = strand, n = 1), "numeric")
  expect_error(dplyr::slice_sample(gr_rle, n = 1, weight_by = seqnames), "numeric")
})
