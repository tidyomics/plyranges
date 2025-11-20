
# Setup test data
setup_test_data <- function() {
  # Create basic GRanges with various column types
  df <- data.frame(
    seqnames = rep(c("chr1", "chr2"), each = 5),
    start = 1:10,
    end = 5:14,
    strand = rep(c("+", "-"), 5),
    score = runif(10),
    count = rpois(10, 5),
    category = factor(rep(c("A", "B"), 5))
  )
  
  rng <- as_granges(df)
  
  # Add S4 columns to mcols for testing
  mcols(rng)$s4_rle <- Rle(c("x", "y", "z"), c(3, 4, 3))
  mcols(rng)$factor_col <- as(factor(rep(c("group1", "group2"), 5)), "Rle")
  
  return(rng)
}

test_that("pull works with core GRanges columns", {
  rng <- setup_test_data()
  
  # Test core genomic columns
  expect_equal(pull(rng, start), start(rng))
  expect_equal(pull(rng, end), end(rng))
  expect_equal(pull(rng, width), width(rng))
  expect_equal(pull(rng, seqnames), seqnames(rng))
  expect_equal(pull(rng, strand), strand(rng))
})

test_that("pull works with metadata columns", {
  rng <- setup_test_data()
  
  # Test regular mcols
  expect_equal(pull(rng, score), mcols(rng)[["score"]])
  expect_equal(pull(rng, count), mcols(rng)[["count"]])
  expect_equal(pull(rng, category), mcols(rng)[["category"]])
})

test_that("pull works with S4 metadata columns", {
  rng <- setup_test_data()
  
  # Test S4 columns in mcols
  expect_equal(pull(rng, s4_rle), mcols(rng)[["s4_rle"]])
  expect_equal(pull(rng, factor_col), mcols(rng)[["factor_col"]])
  
  # Verify they are indeed S4 objects
  expect_true(isS4(mcols(rng)[["s4_rle"]]))
  expect_true(isS4(mcols(rng)[["factor_col"]]))
})

test_that("pull works with positional indexing", {
  rng <- setup_test_data()
  
  # Test positive indexing (metadata columns only, excluding core columns)
  mcol_names <- setdiff(tbl_vars(rng), c("seqnames", "start", "end", "width", "strand"))
  
  expect_equal(pull(rng, 1), start(rng))
  expect_equal(pull(rng, 2), end(rng))
  
  # Test negative indexing (from the right)
  expect_equal(pull(rng, -1), mcols(rng)[[mcol_names[length(mcol_names)]]])
  expect_equal(pull(rng, -2), mcols(rng)[[mcol_names[length(mcol_names) - 1]]])
})

test_that("pull with name parameter works correctly", {
  rng <- setup_test_data()
  
  # Test naming with regular columns
  result <- pull(rng, score, name = category)
  expected <- mcols(rng)[["score"]]
  names(expected) <- mcols(rng)[["category"]]
  expect_equal(result, expected)
  
  # Test naming with core column
  result <- pull(rng, count, name = seqnames)
  expected <- mcols(rng)[["count"]]
  names(expected) <- as.character(seqnames(rng))
  expect_equal(result, expected)
})

test_that("pull with S4 name parameter converts to vector", {
  rng <- setup_test_data()
  
  # Test naming with S4 column - should convert to vector
  result <- pull(rng, score, name = s4_rle)
  expected <- mcols(rng)[["score"]]
  names(expected) <- as.vector(mcols(rng)[["s4_rle"]])
  expect_equal(result, expected)
  
  # Test naming with S4 factor
  result <- pull(rng, count, name = factor_col)
  expected <- mcols(rng)[["count"]]
  names(expected) <- as.vector(mcols(rng)[["factor_col"]])
  expect_equal(result, expected)
})

test_that("pull works with GroupedGenomicRanges", {
  rng <- setup_test_data()
  grouped_rng <- group_by(rng, strand)
  
  # Test that pull on grouped ranges returns same as ungrouped
  expect_equal(pull(grouped_rng, start), start(rng))
  expect_equal(pull(grouped_rng, score), mcols(rng)[["score"]])
  expect_equal(pull(grouped_rng, s4_rle), mcols(rng)[["s4_rle"]])
  
  # Test with names on grouped ranges
  result_grouped <- pull(grouped_rng, score, name = category)
  result_ungrouped <- pull(rng, score, name = category)
  expect_equal(result_grouped, result_ungrouped)
  
  # Test with S4 names on grouped ranges
  result_grouped <- pull(grouped_rng, count, name = s4_rle)
  result_ungrouped <- pull(rng, count, name = s4_rle)
  expect_equal(result_grouped, result_ungrouped)
})

test_that("pull works with GroupedIntegerRanges", {
  # Create IRanges version
  df <- data.frame(
    start = 1:10,
    width = 5,
    score = runif(10),
    group = rep(c("A", "B"), 5)
  )
  
  irng <- as_iranges(df)
  mcols(irng)$s4_rle <- Rle(c("x", "y"), c(5, 5))
  
  grouped_irng <- group_by(irng, group)
  
  # Test basic functionality
  expect_equal(pull(grouped_irng, start), start(irng))
  expect_equal(pull(grouped_irng, score), mcols(irng)[["score"]])
  expect_equal(pull(grouped_irng, s4_rle), mcols(irng)[["s4_rle"]])
  
  # Test with names
  result <- pull(grouped_irng, score, name = group)
  expected <- mcols(irng)[["score"]]
  names(expected) <- mcols(irng)[["group"]]
  expect_equal(result, expected)
})

test_that("pull handles edge cases appropriately", {
  rng <- setup_test_data()
  
  # Test quasiquotation support
  var_name <- "score"
  expect_equal(pull(rng, !!var_name), mcols(rng)[["score"]])
  
  # Test with symbols
  expect_equal(pull(rng, !!rlang::sym("count")), mcols(rng)[["count"]])
})

test_that("pull returns correct data types", {
  rng <- setup_test_data()
  
  # Verify return types match original column types
  expect_identical(class(pull(rng, score)), class(mcols(rng)[["score"]]))
  expect_identical(class(pull(rng, count)), class(mcols(rng)[["count"]]))
  expect_identical(class(pull(rng, category)), class(mcols(rng)[["category"]]))
  expect_identical(class(pull(rng, s4_rle)), class(mcols(rng)[["s4_rle"]]))
  expect_identical(class(pull(rng, start)), class(start(rng)))
  expect_identical(class(pull(rng, seqnames)), class(seqnames(rng)))
})

test_that("pull default argument (var = -1) works correctly", {
  rng <- setup_test_data()
  
  # Test default behavior - should pull last column
  all_vars <- tbl_vars(rng)
  last_var <- all_vars[length(all_vars)]
  
  # Default pull() should equal pull with last column explicitly
  expect_equal(pull(rng), pull(rng, !!last_var))
  
  # Test with grouped ranges
  grouped_rng <- group_by(rng, strand)
  expect_equal(pull(grouped_rng), pull(grouped_rng, !!last_var))
  
  # Test that default pulls from the rightmost position
  expect_equal(pull(rng), pull(rng, -1))
})

test_that("pull positional indexing works correctly", {
  rng <- setup_test_data()
  all_vars <- tbl_vars(rng)
  
  # Test positive indexing from left
  for (i in 1:min(5, length(all_vars))) {
    expect_equal(pull(rng, i), pull(rng, !!all_vars[i]),
                 info = paste("Failed for position", i))
  }
  
  # Test negative indexing from right
  for (i in 1:min(5, length(all_vars))) {
    expected_pos <- length(all_vars) - i + 1
    expect_equal(pull(rng, -i), pull(rng, !!all_vars[expected_pos]),
                 info = paste("Failed for negative position", -i))
  }
  
  # Test with grouped ranges
  grouped_rng <- group_by(rng, strand)
  expect_equal(pull(grouped_rng, 1), pull(rng, 1))
  expect_equal(pull(grouped_rng, -1), pull(rng, -1))
  expect_equal(pull(grouped_rng, 3), pull(rng, 3))
})

test_that("pull positional indexing with different Range types", {
  # Test with IRanges
  df_ir <- data.frame(
    start = 1:5,
    width = 3,
    score1 = runif(5),
    score2 = rpois(5, 3),
    group = letters[1:5]
  )
  irng <- as_iranges(df_ir)
  ir_vars <- tbl_vars(irng)
  
  # Test first and last positions
  expect_equal(pull(irng, 1), pull(irng, !!ir_vars[1]))
  expect_equal(pull(irng, -1), pull(irng, !!ir_vars[length(ir_vars)]))
  expect_equal(pull(irng), pull(irng, -1))  # Default behavior
  
  # Test with minimal IRanges (only core columns)
  minimal_irng <- IRanges(start = 1:3, width = 5)
  minimal_vars <- tbl_vars(minimal_irng)
  expect_equal(pull(minimal_irng, 1), pull(minimal_irng, !!minimal_vars[1]))
  expect_equal(pull(minimal_irng), pull(minimal_irng, -1))
})

test_that("pull default behavior with different column compositions", {
  # Test with ranges that have only core columns
  core_only <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(1:3, width = 5),
    strand = "+"
  )
  core_vars <- tbl_vars(core_only)
  expect_equal(pull(core_only), pull(core_only, !!core_vars[length(core_vars)]))
  
  # Test with ranges that have mixed core + metadata columns
  rng <- setup_test_data()
  all_vars <- tbl_vars(rng)
  expect_equal(pull(rng), pull(rng, !!all_vars[length(all_vars)]))
  
  # Verify the last column is indeed what we expect
  last_col_name <- all_vars[length(all_vars)]
  if (last_col_name %in% names(mcols(rng))) {
    expect_equal(pull(rng), mcols(rng)[[last_col_name]])
  } else {
    # It's a core column
    if (last_col_name == "start") expect_equal(pull(rng), start(rng))
    if (last_col_name == "end") expect_equal(pull(rng), end(rng))
    if (last_col_name == "width") expect_equal(pull(rng), width(rng))
    if (last_col_name == "seqnames") expect_equal(pull(rng), seqnames(rng))
    if (last_col_name == "strand") expect_equal(pull(rng), strand(rng))
  }
})

test_that("pull with empty mcols still accesses core columns", {
  # Create ranges with no metadata columns
  rng <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(start = 1:5, end = 5:9),
    strand = "+"
  )
  
  # Should still be able to pull core columns
  expect_equal(pull(rng, start), start(rng))
  expect_equal(pull(rng, seqnames), seqnames(rng))
  expect_equal(pull(rng, strand), strand(rng))
  
  # Test default behavior with no mcols
  core_vars <- tbl_vars(rng)
  expect_equal(pull(rng), pull(rng, !!core_vars[length(core_vars)]))
  expect_equal(pull(rng), pull(rng, -1))
})

test_that("pull maintains vector attributes", {
  rng <- setup_test_data()
  
  # Test that factor levels are preserved
  original_factor <- mcols(rng)[["category"]]
  pulled_factor <- pull(rng, category)
  expect_equal(levels(original_factor), levels(pulled_factor))
  
  # Test that S4 object structure is preserved
  original_rle <- mcols(rng)[["s4_rle"]]
  pulled_rle <- pull(rng, s4_rle)
  expect_identical(runLength(original_rle), runLength(pulled_rle))
  expect_identical(runValue(original_rle), runValue(pulled_rle))
})