context("unnesting GRanges objects")

test_that("unnesting makes sense", {
  # unnesting empty GRanges returns error
  gr <- GRanges()
  expect_error(unnest(gr), "No list columns to unnest.", fixed = TRUE)
  # no metadata columns returns error
  gr <- as_granges(data.frame(seqnames = c("chr1", "chr2"),
                              start = c(10, 20),
                              width = 5))
  expect_error(unnest(gr), "No list columns to unnest.", fixed = TRUE)

  # by default unnest does cartesian product of list columns
  gr <- as_granges(data.frame(seqnames = "chr1", start = 20:22, width = 1000))
  gr <- mutate(gr,
                 col1 = IntegerList(a = 1, b = c(4,5), c = c(2,3)),
                 col2 = IntegerList(c(1,2), c(3,4), c(5)),
                 score = 1:3)

  correct_gr <- GRanges(seqnames = "chr1",
                        ranges = IRanges(start = c(20,20,21,21,21,21,22,22),
                                         width = 1000),
                        col1 = c(1L,1L,4L,4L,5L,5L,2L,3L),
                        col2 = c(1L,2L,3L,4L,3L,4L,5L,5L),
                        score = c(1L,1L,2L,2L,2L,2L,3L,3L))
  test_gr <- unnest(gr)
  expect_identical(correct_gr, test_gr)
  test_gr <- unnest(gr, col1, col2)
  expect_identical(correct_gr, test_gr)

  # unnesting on non-existent column errors
  expect_error(unnest(gr, gc),
               "Input column(s): gc not found.",
               fixed = TRUE)
  # unnesting on non-list column returns an error
  expect_error(unnest(gr, score),
               "Input column(s): score are not list columns.",
               fixed = TRUE)

  # not specifying .drop returns other list columns
  correct_gr <- GRanges(seqnames = "chr1",
                          ranges = IRanges(start = c(20,21,21,22,22),
                                           width = 1000),
                          col1 = as.integer(c(1,4,5,2,3)),
                          col2 = IntegerList(c(1,2), c(3,4), c(3,4), 5,5),
                          score = as.integer(c(1,2,2,3,3)))
  test_gr <- unnest(gr, col1)
  expect_identical(correct_gr, test_gr)

  correct_gr <- select(correct_gr, -col2)
  test_gr <- unnest(gr, col1, .drop = TRUE)
  expect_identical(correct_gr, test_gr)

  # .id works as expected
  expect_error(unnest(gr, .id = "name"),
               "`.id` does not have same length as number of list columns.",
               fixed = TRUE)
  expect_error(unnest(gr, col2, .id = c("name", "x0")),
               "`.id` does not have same length as number of list columns.",
               fixed = TRUE)

  correct_gr <- GRanges(seqnames = "chr1",
                        ranges = IRanges(start = c(20,20,21,21,21,21,22,22),
                                         width = 1000),
                        col1 = as.integer(c(1,1,4,4,5,5,2,3)),
                        col2 = as.integer(c(1,2,3,4,3,4,5,5)),
                        score = as.integer(c(1,1,2,2,2,2,3,3)),
                        id1 = c(rep("a", 2), rep("b", 4), rep("c", 2)),
                        id2 = as.integer(c(1,1, 2,2,2,2,3,3)))

  test_gr <- unnest(gr, .id = c("id1", "id2"))
  expect_identical(correct_gr, test_gr)

  correct_gr <- GRanges(seqnames = "chr1",
                        ranges = IRanges(start = c(20,21,21,22,22),
                                         width = 1000),
                        col1 = as.integer(c(1,4,5,2,3)),
                        col2 = IntegerList(c(1,2), c(3,4), c(3,4), 5,5),
                        score = as.integer(c(1,2,2,3,3)),
                        id1 = c("a", "b", "b", "c", "c"))
  test_gr <- unnest(gr, col1, .id = "id1")
  expect_identical(correct_gr, test_gr)
  
  # one column test
  gr <- as_granges(data.frame(seqnames = "chr1", start = 20:22, width = 1000))
  gr <- mutate(gr, col1 = IntegerList(a = 1, b = c(4,5), c = c(2,3)))
  correct_gr <- S4Vectors::expand(gr)
  expect_equal(unnest(gr), correct_gr)
  
  correct_gr$id1 <- c("a", "b", "b", "c", "c")
  expect_equal(unnest(gr, .id = "id1"), correct_gr)
  
  # drop emptys
  gr <- as_granges(data.frame(seqnames = "chr1", start = 20:22, width = 1000))
  gr <- mutate(gr, col1 = IntegerList(a = integer(), b = c(4,5), c = c(2,3)))
  correct_gr <- S4Vectors::expand(gr, keepEmptyRows = TRUE)
  expect_equal(unnest(gr, .keep_empty = TRUE), correct_gr)
  correct_gr$id1 <- c("a", "b", "b", "c", "c")
  expect_equal(unnest(gr, .id = "id1", .keep_empty = TRUE), correct_gr)

})
