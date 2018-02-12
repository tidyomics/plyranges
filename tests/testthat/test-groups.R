# test-groups.R
context("Grouped Ranges and grouping")

gr0  <- GRanges()
ir0 <- IRanges()
gr1 <- GRanges(seqnames = "seq1",
               strand = c(rep("+", 4), rep("-", 3), rep("*", 3)),
               ranges = IRanges(start = 1:10, width = 5),
               score = as.integer(c(0,3,4,5,10,14,0, 0, 4,8)))

test_that("grouping works on empty ranges", {
  expect_s4_class(gr0 %>% group_by(seqnames), "GRangesGrouped")
  expect_s4_class(ir0 %>% group_by(start), "IRangesGrouped")
})

test_that("throws invalid object if group not found", {
  expect_error(gr0 %>% group_by(st))
  expect_error(ir0 %>% group_by(st))
})

test_that("grouping allows character names", {
  expect_s4_class(gr0 %>% group_by("seqnames"), "GRangesGrouped")
  expect_s4_class(gr1 %>% group_by("score"), "GRangesGrouped")
  expect_identical(gr1 %>% group_by(score), gr1 %>% group_by("score"))
})

test_that("group_vars works as expected", {
  gr_by_strand_score <- gr1 %>% group_by(strand, score)
  expect_equal(group_vars(gr_by_strand_score),
                   c("strand", "score"))
  ir_by_start <- ir0 %>% group_by(start)
  expect_equal(group_vars(ir_by_start), "start")
})

test_that("ungrouping works as expected", {
  gr_by_strand_score <- gr1 %>% group_by(strand, score)
  gr_by_score <- gr1 %>% group_by(score)
  expect_identical(ungroup(gr_by_strand_score), gr1)
  expect_equivalent(gr_by_strand_score %>% ungroup(strand), gr_by_score)
})

test_that("group by matches HelloRanges", {
  skip_on_travis()
  oldwd <- getwd()
  setwd(system.file("unitTests", "data", "groupby", package="HelloRanges"))

  a <- read_bed("values3.header.bed")
  exp <- S4Vectors::aggregate(unstrand(a), score.sum = sum(score)) %>%
    select(-grouping)

  result <- a %>%
    group_by(seqnames, start, end) %>%
    summarise(score.sum = sum(score)) %>%
    as_granges() %>%
    sort()

  expect_identical(exp, result)
  setwd(oldwd)
})
