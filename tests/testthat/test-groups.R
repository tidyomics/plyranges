# test-groups.R
context("Grouped Ranges and grouping")

gr0  <- GRanges()
ir0 <- IRanges()
gr1 <- GRanges(seqnames = "seq1",
               strand = c(rep("+", 4), rep("-", 3), rep("*", 3)),
               ranges = IRanges(start = 1:10, width = 5),
               score = as.integer(c(0,3,4,5,10,14,0, 0, 4,8)))

test_that("grouping works on empty ranges", {
  expect_s4_class(gr0 %>% group_by(seqnames), "GroupedGenomicRanges")
  expect_s4_class(ir0 %>% group_by(start), "GroupedIntegerRanges")
})

test_that("throws invalid object if group not found", {
  expect_error(gr0 %>% group_by(st))
  expect_error(ir0 %>% group_by(st))
})

test_that("grouping doesn't allow character names", {
  expect_error(gr0 %>% group_by("seqnames"))
  expect_error(ir0 %>% group_by("start"))
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

test_that("grouping after filter or mutate leaves groups in tact", {
            by_strand <- gr1 %>% group_by(strand) %>% mutate(new_score = mean(score))
            expect_s4_class(by_strand, "GroupedGenomicRanges")
            expect_equal("strand", group_vars(by_strand))
            by_strand_filter <- gr1 %>% group_by(score) %>% filter(n() >= 2)
            expect_s4_class(by_strand_filter, "GroupedGenomicRanges")
            expect_equal("score", group_vars(by_strand_filter))
})

test_that("group by matches HelloRanges", {
  oldwd <- getwd()
  setwd(system.file("unitTests", "data", "groupby", package="HelloRanges"))

  a <- read_bed("values3.header.bed")
  exp <- GRanges(seqnames = c(rep("chr1", 5), rep("chr3",4)),
                 ranges = IRanges(start = c(1, 11, 12, 21, 121,1,11,21,121),
                                  end = c(10, 20, 21, 30, 130, 10, 20, 30, 130)),
                 strand = "*",
                 score.sum = c(10,5,5,45,1,1,2,3,8))

  result <- a %>%
    group_by(seqnames, start, end) %>%
    summarise(score.sum = sum(score)) %>%
    as_granges() %>%
    sort()

  expect_identical(exp, result)
  setwd(oldwd)
})

test_that("ungroup/groups/group_vars on ordinary Ranges",{
  expect_identical(ungroup(gr1), gr1)
  expect_identical(ungroup(ranges(gr1)), ranges(gr1))
  expect_null(groups(gr1))
  expect_null(groups(ranges(gr1)))
  expect_identical(group_vars(gr1), character(0))
  expect_identical(group_vars(ranges(gr1)), character(0))
})

# issue #42
test_that("group_by does not allow unknown columns", {
  gr <- GRanges("chr1", IRanges(start = 1:6, width = 10))
  foo <- 1:6
  expect_error(group_by(gr, foo), "Column `foo` is unknown")
  expect_error(group_by(gr, strand, foo), "Column `foo` is unknown")
})



