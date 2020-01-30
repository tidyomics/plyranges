context("disjoin ranges")

test_that("matches IRanges/GRanges tests", {
  x <- IRanges()
  expect_identical(x, disjoin_ranges(x))
  ir <- IRanges(c(1, 21, 10, 1, 15, 5, 20, 20),
                c(6, 20, 9, 3, 14, 11, 20, 19))
  correct_ir <- IRanges(c(1, 4, 5, 7, 10, 20), c(3, 4, 6, 9, 11, 20))
  test_ir <- disjoin_ranges(ir)
  expect_identical(test_ir, correct_ir)

  # check revmap
  test_ir <- ir %>%
    mutate(i = 1:n()) %>%
    disjoin_ranges(revmap = IRanges::IntegerList(i))

  mcols(correct_ir)$revmap <- IRanges::IntegerList(c(1, 4), 1, c(1, 6), 6, 6, 7)
  expect_identical(test_ir, correct_ir)

  # -- granges
  gr <- GRanges(Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
          IRanges(1:10, width=10:1, names=head(letters, 10)),
          Rle(c("-", "+", "*", "+", "-"), c(1, 2, 2, 3, 2)),
          score=1:10, GC=seq(1, 0, length=10),
          seqinfo=Seqinfo(paste("chr", 1:3, sep="")))

  correct_gr <- GRanges(Rle(c("chr1", "chr2", "chr3"), c(3, 3, 4)),
                        IRanges(start=c(6, 1, 5, 2, 3, 4, 7, 8, 9, 10),
                                end=c(10, 10, 10, 2, 10, 10, 7, 10, 9, 10)),
                        c("+", "-", "*", "+", "+", "*", "+", "+", "-", "-"))
  # matches directed
  expect_identical(disjoin_ranges_directed(gr), correct_gr)
  # this is the same as disjoin unstranded on correct_gr
  expect_identical(disjoin_ranges(gr), disjoin_ranges(correct_gr))

  gr <- GRanges(Rle(c("chr1", "chr3"), c(2, 2)),
                IRanges(c(8, 6, 8, 6), c(11, 15, 11, 15),
                        names=c("k", "l", "m", "n")),
                c("-", "-", "+", "*"),
                score=11:14, GC=c(.2, .3, .3, .1))

  correct_gr <- GRanges(Rle(c("chr1", "chr3"), c(3, 2)),
                        IRanges(c(6, 8, 12, 8, 6), c(7, 11, 15, 11, 15)),
                        Rle(c("-", "+", "*"), c(3, 1, 1)))
  mcols(correct_gr)$revmap <- IRanges::IntegerList(2, 1:2, 2, 3, 4)
  expect_identical(gr %>%
                     mutate(i = 1:n()) %>%
                     disjoin_ranges_directed(revmap = IRanges::IntegerList(i)),
                   correct_gr)

  # grouping works as expected
  grl <- GRangesList( GRanges(Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
                              IRanges(1:10, width=10:1),
                              Rle(c("-", "+", "*", "+", "-"), c(1, 2, 2, 3, 2)),
                              score=1:10, GC=seq(1, 0, length=10),
                              seqinfo=Seqinfo(paste("chr", 1:3, sep=""))),
                      GRanges("1", IRanges(1, 10), score=21, GC=.21),
                      GRanges(),
                      GRanges(Rle(c("chr1", "chr3"), c(2, 2)),
                              IRanges(c(8, 6, 8, 6), c(11, 15, 11, 15)),
                              strand(c("-", "-","+","*")),
                              score=41:44, GC=c(.41, .42, .43, .44)))

  gr_by_group <- stack(grl, "name") %>% group_by(name)

  target <- stack(disjoin(grl, ignore.strand = TRUE), "name")
  current <- disjoin_ranges(gr_by_group) %>% 
    mutate(name = Rle(name))
  expect_identical(target, current)
})



test_that("matches HelloRanges multinter", {
  oldwd <- getwd()
  setwd(system.file("unitTests", "data", "multiinter", package="HelloRanges"))
  bed_files <- list.files(pattern = ".bed$")

  correct_gr <- GRanges("chr1",
                        IRanges(c(7, 9, 13, 16, 21, 23, 31, 33),
                                c(8, 12, 15, 20, 22, 30, 32, 34)),
                        i=IRanges::FactorList(1, c(1,3), 1:3, 1:2, 2, 1:2, 2, 3))

  gr_l <- S4Vectors::List(lapply(bed_files, function(x) {
    mutate(read_bed(x), grp = sub(".bed$", "", basename(x)))
    }))

  gr_by_group_r <- unlist(gr_l) %>%
    mutate(grp = factor(grp, levels = c("a", "b", "c"))) %>%
    group_by(grp) %>%
    reduce_ranges()
  test_gr <- gr_by_group_r %>%
    mutate(i = factor(as.integer(grp))) %>%
    disjoin_ranges(i = IRanges::FactorList(i))

  expect_identical(correct_gr, test_gr)

  # with names in place of integer
  mcols(correct_gr)$i <- extractList(factor(c("a", "b", "c")),
                              IRanges::IntegerList(mcols(correct_gr)$i))
  test_gr <- gr_by_group_r %>%
    disjoin_ranges(i = IRanges::FactorList(grp))
  expect_identical(correct_gr, test_gr)

  setwd(oldwd)

})
