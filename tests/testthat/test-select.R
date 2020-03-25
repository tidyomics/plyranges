context("Column wise selection")

gr0 <- GRanges()
ir0 <- IRanges()

mdata <- DataFrame(score = as.integer(c(0,3,4,5,10,14,0, 0, 4,8)),
                   gc = runif(10),
                   counts = rpois(10, 2))
ir1 <- IRanges(start = 1:10, width = 5)
mcols(ir1) <- mdata
gr1 <- GRanges(seqnames = "seq1",
               strand = c(rep("+", 4), rep("-", 3), rep("*", 3)),
               ranges = ir1)


test_that("no selection returns the same object", {
  expect_identical(select(ir0), ir0)
  expect_identical(select(gr0), gr0)
  expect_identical(select(ir1), ir1)
  expect_identical(select(gr1), gr1)
})

test_that("dropping core parts throws an error...", {
  expect_error(select(ir0, start))
  expect_error(select(gr0, start))
})

test_that("...unless allowing to drop ranges", {
  res <- select(ir0, start, .drop_ranges = TRUE)
  
  expect_s4_class(res, "DataFrame")
  expect_equal(ncol(res), 1)
  expect_equal(length(res$start), 0)
  
  res <- select(gr0, start, .drop_ranges = TRUE)
  expect_s4_class(res, "DataFrame")
  expect_equal(ncol(res), 1)
  expect_equal(length(res$start), 0)
  
  expect_identical(select(ir1, score, gc, counts, .drop_ranges = TRUE),
                   mcols(ir1))
  expect_identical(select(gr1, score, gc, counts, .drop_ranges = TRUE),
                   mcols(gr1))

})

test_that("reordering/dropping works as expected", {
  gr2 <- select(gr1, counts, gc, score)
  expect_equal(names(mcols(gr2)), c("counts", "gc", "score"))
  gr3 <- select(gr1, -1)
  gr4 <- select(gr1, -score)
  expect_identical(gr3, gr4)
  ir2 <- select(ir1, -counts)
  expect_equal(names(mcols(ir2)), c("score", "gc"))
  ir3 <- select(ir1, 1:2)
  expect_identical(ir2, ir3)
})

test_that("dropping everything sets mcols slot to empty", {
  gr5 <- select(gr1, -score:-counts)
  mcols(gr1) <- NULL
  expect_identical(gr5, gr1)
  ir4 <- select(ir1, -score:-counts)
  expect_equivalent(ir4, IRanges(ir1))
})

