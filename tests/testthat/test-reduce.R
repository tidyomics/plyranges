context("reduce_ranges method")

test_that("matches IRanges/GenomicRanges", {
  x <- IRanges()
  expect_identical(x, reduce_ranges(x))

  x <- IRanges(c(1:4, 10:11, 11), width=c(0,1,1,0,0,0,1))
  mcols(x) <- DataFrame(mapping = paste0("a", seq_along(x)))
  target <- IRanges(c(1:2, 10:11), width=c(0,2,0,1))
  mcols(target) <- DataFrame(mapping=c("a1","a2,a3,a4", "a5", "a6,a7"))
  expect_identical(reduce_ranges(x, mapping = paste0(mapping, collapse = ",")),
                   target)

  # drop.empty.ranges is just a filter
  current <- x %>%
    filter(width > 0) %>%
    reduce_ranges(mapping = paste0(mapping, collapse = ","))
  target <- reduce(x, drop.empty.ranges=TRUE)
  mcols(target) <- DataFrame(mapping = c("a2,a3", "a7"))

  expect_identical(current, target)

  # -- GRanges
  gr <- GRanges(Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
          IRanges(1:10, width=10:1, names=head(letters, 10)),
          Rle(c("-", "+", "*", "+", "-"), c(1, 2, 2, 3, 2)),
          name = paste0("a", 1:10))

  target <- GRanges(Rle(c("chr1", "chr2", "chr3"), c(3, 2, 2)),
                    IRanges(start=c(6, 1, 5, 2, 4, 7, 9), end=10),
                    c("+", "-", "*", "+", "*", "+", "-"))
  expect_identical(reduce_ranges_directed(gr), target)

  mcols(target)$mapping <- c("a6", "a1", "a5", "a2,a3", "a4", "a7,a8", "a9,a10")

  expect_identical(reduce_ranges_directed(gr, mapping = paste0(name, collapse = ",")),
                   target)
  target <- GRanges(Rle(c("chr1", "chr2", "chr3"), c(1,1,1)),
                    IRanges(start = c(1,2,7), end = c(10,10,10)))
  expect_identical(reduce_ranges(gr), target)

  mcols(target)$mapping <- c("a1,a5,a6", "a2,a3,a4", "a7,a8,a9,a10")

  expect_identical(reduce_ranges(gr, mapping = paste0(name, collapse = ",")),
                   target)


})


test_that("non-standard evaluation works as expected",{
  oldwd <- getwd()
  setwd(system.file("unitTests", "data", "merge", package="HelloRanges"))

  a <- read_bed("a.bed")

  exp <- reduce(a, ignore.strand = TRUE)

  expect_identical(exp, reduce_ranges(a))

  exp <- reduce(a, with.revmap=TRUE, ignore.strand = TRUE)
  mcols(exp) <- aggregate(a, mcols(exp)$revmap,
                          seqnames.count = lengths(seqnames))
  exp <- exp %>% select(-grouping)

  expect_identical(exp, a %>% reduce_ranges(seqnames.count = length(seqnames)))

  mcols(a)$name <- paste0("a", 1:4)
  exp <- reduce(a, with.revmap=TRUE)
  mcols(exp) <- aggregate(a, mcols(exp)$revmap,
                          name.collapse = unstrsplit(name, ","))
  exp <- exp %>% select(-grouping)

  expect_identical(exp, a %>%
                     reduce_ranges(name.collapse = paste0(name, collapse = ","))
                   )
  a <- read_bed("a.full.bed")
  exp <- reduce(a, with.revmap=TRUE, ignore.strand=TRUE)
  mcols(exp) <- aggregate(a, mcols(exp)$revmap,
                          name.collapse = unstrsplit(name, ","),
                          score.sum = sum(score))
  exp <- exp %>% select(-grouping)

  expect_identical(exp, reduce_ranges(a,
                                      name.collapse = paste0(name, collapse = ","),
                                      score.sum = sum(score)))

  exp <- reduce(a, with.revmap=TRUE, ignore.strand=TRUE)
  mcols(exp) <- aggregate(a, mcols(exp)$revmap,
                          score.count = lengths(score),
                          score.sum = sum(score))
  exp <- exp %>% select(-grouping)

  expect_identical(exp, reduce_ranges(a, score.count = length(score),
                                      score.sum = sum(score)))
  setwd(oldwd)
})
