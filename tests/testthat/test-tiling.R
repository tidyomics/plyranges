context("tile_ranges and slide_ranges")

gr <- GRanges(
  seqnames = "1",
  ranges   = IRanges(start = c(1, 101), width = c(60, 80)),
  strand   = c("+", "-")
)

test_that("tile_ranges returns partition column and correct ranges", {
  res <- tile_ranges(gr, width = 20)
  expect_true(!is.null(mcols(res)[["partition"]]))
})

test_that("slide_ranges returns partition column and correct ranges", {
  res <- slide_ranges(gr, width = 20, step = 10)
  expect_true(!is.null(mcols(res)[["partition"]]))
})

test_that("tile_ranges_directed reverses negative strand tiles", {
  fwd <- tile_ranges(gr, width = 20)
  dir <- tile_ranges_directed(gr, width = 20)
  neg_fwd <- fwd[mcols(fwd)$partition == 2L]
  neg_dir <- dir[mcols(dir)$partition == 2L]
  expect_identical(neg_fwd, rev(neg_dir))
})

test_that("slide_ranges_directed reverses negative strand windows", {
  fwd <- slide_ranges(gr, width = 20, step = 10)
  dir <- slide_ranges_directed(gr, width = 20, step = 10)
  neg_fwd <- fwd[mcols(fwd)$partition == 2L]
  neg_dir <- dir[mcols(dir)$partition == 2L]
  expect_identical(neg_fwd, rev(neg_dir))
})
