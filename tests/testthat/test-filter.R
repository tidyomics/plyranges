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
  expect_identical(gr0 %>% filter(score > 0, strand == "*"),
                   gr0[score(gr0) > 0 & strand(gr0) == "*"])
})

test_that("checking calls to n", {
  expect_identical(gr0 %>% filter(n() == 10), gr0)
  expect_identical(gr0 %>% filter(n() == 4), GRanges(score = integer(0),
                                                     seqinfo = seqinfo(gr0)))
  expect_error(gr0 %>% filter(1:n()))
})

gr_gfilter <- c(gr0[strand(gr0) == "+" & score(gr0) > 3],
                gr0[strand(gr0) == "-" & score(gr0) > 8],
                gr0[strand(gr0) == "*" & score(gr0) > 4])
test_that("grouped filter checks", {
  expect_identical(gr0 %>%
                     group_by(strand) %>%
                     filter(n() == 4) %>%
                     ungroup(),
                   gr0[1:4])
  expect_identical(gr0 %>% group_by(strand) %>%
                     filter(score > mean(score)) %>%
                     ungroup(),
                   gr_gfilter)
})

# see https://github.com/sa-lee/plyranges/issues/36
test_that("using `::` works", {
  expect_identical(gr0 %>% plyranges::group_by(strand) %>%
                     plyranges::filter(score > mean(score)) %>%
                     plyranges::ungroup(),
                   gr_gfilter)
  gr <- data.frame(seqnames = c("chr1", "chr2", "chr2", "chr1", "chr2"),
                   start = 1:5,
                   width = 5) %>%
    as_granges()
  expect_identical(plyranges::filter(gr, seqnames %in% c("chr1", "chr2")),
                   gr)
  expect_identical(gr %>% 
                     plyranges::group_by(seqnames) %>% 
                     plyranges::filter(seqnames %in% c("chr1", "chr2")) %>%
                     plyranges::ungroup(),
                   gr)
})
