context("Anchored ranges")

test_that("IRangesAnchored", {
  ir <- IRanges()
  astart <- anchor_start(ir)
  # check class and accessors
  expect_s4_class(astart, "IRangesAnchored")
  expect_equal(anchor(astart), "start")
  # adding a new anchored overrides existing one
  expect_equal(anchor(anchor_center(astart)), "center")
  # can't anchor by strand for IRanges
  expect_error(anchor_5p(astart))
  expect_error(new("IRangesAnchored", ir, anchor = "3p"))
  # stretching empty ranges returns identical range if anchored
  expect_identical(ir, stretch(astart, 5L))

  # anchoring coordinates leaves respective positions fixed
  ir <- IRanges(start = c(1, 20, 25, 25, 33),
                width = c(19, 5, 0, 8, 5))
  correct_start <- start(ir)
  test_start <- ir %>% anchor_start() %>% stretch(10L) %>% start()
  expect_equal(test_start, correct_start)
  correct_end <- end(ir)
  test_end <- ir %>% anchor_end() %>% stretch(10L) %>% end()
  expect_equal(test_end, correct_end)
  correct_center <- resize(ir, fix = "center", width = 11L)
  test_center <- set_width(anchor_center(ir), 11L)
  expect_identical(correct_center, test_center)
})

test_that("GRangesAnchored", {

  gr <- GRanges()
  # check class and acessors
  expect_s4_class(anchor_5p(gr), "GRangesAnchored")
  by_3p <- anchor_3p(gr)
  expect_equal(anchor(by_3p), "3p")
  # unanchored returns null
  expect_null(anchor(gr))

  # anchoring coordinates works as expected
  gr <- GRanges(
    seqnames=Rle(paste("chr", c(1, 2, 1, 3), sep=""), c(1, 3, 2, 4)),
    ranges=IRanges(1:10, width=10:1, names=letters[1:10]),
    strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
    score=1:10,
    GC=seq(1, 0, length=10)
  )

  # modifying width anchors by start
  correct_gr <- gr
  width(correct_gr) <- 5L
  expect_identical(set_width(anchor_start(gr), 5L), correct_gr)

  # anchoring by 5p fixes start for negative strand
  expect_identical(set_width(anchor_5p(gr), 5L),
                   resize(gr, fix = "start", width = 5L))
  # anchoring by 3p fixes end for neg strand
  expect_identical(set_width(anchor_3p(gr), 5L),
                   resize(gr, fix = "end", width = 5L))
})
