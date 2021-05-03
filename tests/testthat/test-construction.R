# test-construction.R
context("Ranges construction")

test_that("invalid data.frame call throws an error", {
  expect_error(as_iranges(data.frame()))
  expect_error(as_iranges(data.frame(start = 1)))
})

test_that("valid data.frame returns expected class", {
  df <- data.frame(start = 1:3, width = 2:4)
  expect_s4_class(as_iranges(df), "IRanges")
  df2 <-df
  df2$strand <- "*"
  expect_s4_class(as_iranges(df2), "IRanges")

  # but there's not enough info for GRanges
  expect_error(as_granges(df2))

  df3 <- df2
  df3$seqnames <- "chr1"
  expect_s4_class(as_granges(df3), "GRanges")

  df4 <- df
  df4$seqnames <- "chr1"
  expect_s4_class(as_granges(df4), "GRanges")

})

test_that("non-standard evaluation works as expected", {
  df <- data.frame(st= 1:3, width = 2:4)
  expect_s4_class(as_iranges(df, start = st), "IRanges")
  # adding new column throws an error
  expect_error(as_iranges(df, start = st, gc = 1))
  expect_s4_class(as_granges(df, start = st, seqnames = "chr1"), "GRanges")
})

test_that("out of bounds Ranges throws error", {
  df <- data.frame(start = 1:3, width = 2:4)
  expect_error(as_iranges(df, end = 4))
  df2 <- df
  df2$seqnames <- "1"
  df2$strand <- "negative"
  expect_error(as_granges(df2))
})

test_that("DataFrame input", {
  # tests issue 62 https://github.com/sa-lee/plyranges/issues/62
  df <- DataFrame(st = 1:3, width = 2:4,
                  grps = IRanges::FactorList("a", c("b", "c"), "d", compress = FALSE))

  ir <- as_iranges(df, start = st)
  expect_identical(mcols(ir)$grps, df$grps)

  df$strand <- "+"
  df$chr <- "chr1"
  gr <- as_granges(df, start = st, seqnames = chr)
  expect_identical(mcols(gr)$grps, df$grps)
})

test_that("Working with names", {
  ir <- IRanges::IRanges(start = 1:3, width = 4, names = c("a", "b", "c"))
  expect_identical(unname(ir), remove_names(ir))
  ir_noname <- names_to_column(ir)
  expect_null(names(ir_noname))
  expect_identical(names(ir), mcols(ir_noname)[["name"]])
  ir_id <- id_to_column(ir_noname)
  expect_identical(seq_len(length(ir)), mcols(ir_id)[["id"]])


  gr <- GenomicRanges::GRanges(seqnames = "chr1", ranges = ir)
  expect_identical(unname(gr), remove_names(gr))
  gr_noname <- names_to_column(gr)
  expect_null(names(gr_noname))
  expect_identical(names(gr), mcols(gr_noname)[["name"]])
  gr_id <- id_to_column(gr_noname)
  expect_identical(seq_len(length(gr)), mcols(gr_id)[["id"]])
  expect_identical(mcols(gr_id),
                   S4Vectors::DataFrame(id = seq_len(length(gr)),
                                        name = names(gr)))
})
