# test-construction.R
context("Ranges consturction")

test_that("invalid data.frame call throws an error", {
  expect_error(Ranges(data.frame()))
  expect_error(Ranges(data.frame(start = 1)))
})

test_that("valid data.frame returns expected class", {
  df <- data.frame(start = 1:3, width = 2:4)
  expect_s4_class(Ranges(df), "IRanges")
  df2 <-df
  df2$strand <- "*"
  expect_s4_class(Ranges(df2), "IRanges")

  df3 <- df2
  df3$seqnames <- "chr1"
  expect_s4_class(Ranges(df3), "GRanges")

  df4 <- df
  df4$seqnames <- "chr1"
  expect_s4_class(Ranges(df4), "GRanges")

})

test_that("non-standard evaluation works as expected", {
  df <- data.frame(st= 1:3, width = 2:4)
  expect_s4_class(Ranges(df, start = st), "IRanges")
  expect_s4_class(Ranges(df, start = st, strand = "*"), "IRanges")
  expect_s4_class(Ranges(df, start = st, seqnames = "chr1"), "GRanges")
})

test_that("out of bounds Ranges throws error", {
  df <- data.frame(start = 1:3, width = 2:4)
  expect_error(Ranges(df, end = 4))
  df2 <- df
  df2$seqnames <- "1"
  df2$strand <- "negative"
  expect_error(Ranges(df2))
})
