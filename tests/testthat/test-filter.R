# test-filter.R
context("filter methods")

gr0 <- GRanges(seqnames = "seq1",
               strand = c(rep("+", 4), rep("-", 3), rep("*", 3)),
               ranges = IRanges(start = 1:10, width = 5),
               score = as.integer(c(0,3,4,5,10,14,0, 0, 4,8)))

test_that("simple filter checks", {

  expect_identical(gr0 %>% filter(strand == "+"), gr0[strand(gr0) == "+"])
  expect_identical(gr0 %>% filter(score > 0), gr0[score(gr0) > 0])
  expect_identical(gr0 %>% filter(score > 0 & strand == "*"),
                   gr0[score(gr0) > 0 & strand(gr0) == "*"])
})

gr_gfilter <- c(gr0[strand(gr0) == "+" & score(gr0) > 3],
                gr0[strand(gr0) == "-" & score(gr0) > 8],
                gr0[strand(gr0) == "*" & score(gr0) > 4])
test_that("grouped filter checks", {
  expect_identical(gr0 %>% group_by(score) %>% filter(score > mean(score)),
                   gr_gfilter)
})
