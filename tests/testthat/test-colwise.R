context("colwise mutate/summarise")

mdata <- DataFrame(score = as.integer(c(0,3,4,5,10,14,0, 0, 4,8)),
                   grp = rep(c("A", "B"), 5))

ir1 <- IRanges(start = 1:10, width = 5)
mcols(ir1) <- mdata
gr1 <- GRanges(seqnames = "seq1",
               strand = c(rep("+", 4), rep("-", 3), rep("*", 3)),
               ranges = ir1)

test_that("summarise evaluates correctly", {
  df <- data.frame(mean = 4.8,  n = 10L)
  expect_equal(summarise(ir1, mean = mean(score), n = n()), df)
  expect_equal(summarise(gr1, mean = mean(score), n = n()), df)

  gdf <- data.frame(grp = c("A", "B"), score = c(3.6, 6),stringsAsFactors = FALSE)
  expect_equal(ir1 %>% group_by(grp) %>% summarise(score = mean(score)), gdf)
  expect_equal(gr1 %>% group_by(grp) %>% summarise(score = mean(score)), gdf)

  gdf <- data.frame(grp = c("A", "B"), n = c(5L,5L), stringsAsFactors = FALSE)
  expect_equal(ir1 %>% group_by(grp) %>% summarise(n = n()), gdf)
  expect_equal(gr1 %>% group_by(grp) %>% summarise(n = n()), gdf)
})

test_that("mutate allows out of scope columns", {
  gr2 <- gr1
  mcols(gr2)$score2 <- score(gr2) + 1L
  mcols(gr2)$score3 <- mcols(gr2)$score2 * 2L
  expect_identical(gr1 %>%
                     mutate(score2 = score + 1L,
                            score3 = score2*2L),
                   gr2)

  ir2 <- ir1
  mcols(ir2)$score2 <- mcols(ir2)$score*4L
  mcols(ir2)$score3 <- mcols(ir2)$score2 + 3L
  expect_identical(ir1 %>%
                     mutate(score2 = score * 4L,
                            score3 = score2 + 3L), ir2)

})

test_that("mutating by groups", {
  gr2 <- gr1
  mcols(gr2)$gt_grp_score <- (score(gr1) > 3.6 & mcols(gr1)$grp == "A") |
    (score(gr1) > 6 & mcols(gr1)$grp == "B")
  gr2 <- gr2[order(gr2$grp)]
  # this is still buggy in the sense it returns the results ordered
  expect_identical(gr1 %>%
                     group_by(grp) %>%
                     mutate(gt_grp_score = score > mean(score)) %>%
                     ungroup(),
                   gr2)
  ir2 <- ir1
  mcols(ir2)$lt_grp_score <- (mcols(ir1)$score < 3.6 & mcols(ir1)$grp == "A") |
    (mcols(ir1)$score < 6 & mcols(ir1)$grp == "B")
  ir2 <- ir2[order(mcols(ir2)$grp)]
  expect_identical(ir1 %>%
                     group_by(grp) %>%
                     mutate(lt_grp_score = score < mean(score)) %>%
                     ungroup(),
                   ir2)

})
