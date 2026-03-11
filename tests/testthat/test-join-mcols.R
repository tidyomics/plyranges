test_that("join_mcols works with different types of y", {

  x <- as_granges(data.frame(seqnames = "chr1", start = 1:5, width=10, id=c(1:4,4), foo=letters[1:5]))
  y <- DataFrame(id=c(1:2,2,4), bar=letters[6:9])
  res1 <- x |> join_mcols_left(y, by="id")

  ydf <- data.frame(id=c(1:2,2,4), bar=letters[6:9])
  res2 <- x |> join_mcols_left(ydf, by="id")

  ytib <- tibble(id=c(1:2,2,4), bar=letters[6:9])
  res3 <- x |> join_mcols_left(ytib, by="id")

  expect_equal(res1,res2)
  expect_equal(res1,res3)

  # check inner vs left
  res_inner <- x |> join_mcols_inner(y, by="id")

  expect_equal(res_inner, res1[!is.na(res1$bar)])

})
