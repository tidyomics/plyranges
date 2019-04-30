context("Combining Ranges objects")

gr <- as_granges(data.frame(start = 10:15,
                            width = 5,
                            seqnames = "seq1"))

gr2 <- as_granges(data.frame(start = 11:14,
                            width = 1:4,
                            seqnames = "seq1"))

gr3 <- gr
names(gr3) <- letters[1:5]

test_that("simple bind is same as `c`", {
  expect_identical(bind_ranges(gr), gr)
  expect_identical(bind_ranges(gr, gr2), c(gr, gr2))
  expect_identical(bind_ranges(gr, list(gr, gr2), gr2), c(gr,gr,gr2,gr2))
})

test_that("setting .id column", {
  gr4  <- bind_ranges(list(a = gr, b = gr2), .id = "origin")
  expect_s4_class(gr4$origin, "Rle")
  expect_equal(runLength(gr4$origin), c(length(gr), length(gr2)))
  expect_identical(gr4, bind_ranges(a = gr, b = gr2, .id = "origin"))
  expect_identical(gr4, bind_ranges(a = gr, list(b = gr2), .id = "origin"))
})

test_that("bind preserves names", {
  expect_null(names(bind_ranges(gr, gr2)))
  expect_equal(names(bind_ranges(gr3)), names(gr3))
  gr5 <- bind_ranges(gr3, gr2)
  expect_equal(names(gr5), c(names(gr3), rep("", length(gr2))))
})

