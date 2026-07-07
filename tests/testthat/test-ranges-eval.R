# test-ranges-eval.R
context("NSE with mcols columns named like functions")

gr_gc <- GRanges(
  seqnames = "chr1",
  ranges = IRanges(start = 1:6, width = 10),
  strand = c("+", "+", "-", "-", "*", "*"),
  gc = c(0.1, 0.9, 0.2, 0.8, 0.25, 0.7)
)

test_that("filter sees the mcols column, not the function of the same name", {
  # ground truth: the values in the column, evaluated directly
  truth <- gr_gc[mcols(gr_gc)$gc < 0.3]
  expect_identical(gr_gc %>% filter(gc < 0.3), truth)
})

test_that("mutate reads a function-named mcols column", {
  out <- gr_gc %>% mutate(gc2 = gc * 10)
  expect_equal(mcols(out)$gc2, mcols(gr_gc)$gc * 10)
})

test_that("summarise reads a function-named mcols column", {
  out <- gr_gc %>% summarise(mgc = mean(gc))
  expect_equal(out$mgc, mean(mcols(gr_gc)$gc))
})

test_that("select can pick a function-named mcols column", {
  # select_rng uses parallelVectorNames() to identify the core slots; the
  # S4Vectors #140 made it include mcols names, so selecting a
  # metadata column errored with "Cannot select/rename the following columns".
  expect_identical(mcols(gr_gc %>% select(gc))$gc, mcols(gr_gc)$gc)
  expect_identical(names(mcols(gr_gc %>% select(gc))), "gc")
  # negative selection and .drop_ranges paths
  expect_identical(names(mcols(gr_gc %>% select(-gc))), character(0))
})

test_that("fluentGenomics names_to_column works with a 'baseMean' column", {
  # Mirrors fluentGenomics.Rmd [results-GRanges]: a DESeq2-style GRanges with a
  # 'baseMean' metadata column, promoted to a column via names_to_column().
  # 'baseMean' has no function of that name, so makeFixedColumnEnv() errored.
  gr <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(start = c(10, 20, 30), width = 5),
    baseMean = c(100, 5, 50)
  )
  names(gr) <- c("geneA", "geneB", "geneC")

  out <- names_to_column(gr, "gene_id")
  expect_true("gene_id" %in% names(mcols(out)))
  expect_identical(mcols(out)$gene_id, c("geneA", "geneB", "geneC"))
  expect_equal(mcols(out)$baseMean, c(100, 5, 50))

  # and the direct mutate/filter the vignette relies on
  expect_equal(
    mcols(gr %>% mutate(lb = log2(baseMean)))$lb,
    log2(c(100, 5, 50))
  )
  expect_identical(gr %>% filter(baseMean > 10), gr[c(1, 3)])
})

test_that("grouped NSE reads a function-named mcols column", {
  truth <- c(
    gr_gc[strand(gr_gc) == "+" & mcols(gr_gc)$gc < 0.5],
    gr_gc[strand(gr_gc) == "-" & mcols(gr_gc)$gc < 0.5],
    gr_gc[strand(gr_gc) == "*" & mcols(gr_gc)$gc < 0.5]
  )
  out <- gr_gc %>%
    group_by(strand) %>%
    filter(gc < 0.5) %>%
    ungroup()
  expect_identical(out, truth)

  gout <- gr_gc %>%
    group_by(strand) %>%
    summarise(mgc = mean(gc))
  # one row per strand level present, means computed from the column
  expect_equal(
    sort(gout$mgc),
    sort(sapply(split(mcols(gr_gc)$gc, strand(gr_gc)), mean)),
    check.names = FALSE
  )
})