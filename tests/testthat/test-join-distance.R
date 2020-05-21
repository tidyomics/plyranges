test_nearest_join_distance <- function(query, subject, join_fun, range_fun = granges){
  join <- join_fun(query, subject)
  
  join_dist <- join_fun(query, subject, distance = TRUE)
  join_test <- join_fun(query, subject, distance = "test")
  
  expect_equal(range_fun(join), range_fun(join_dist))
  expect_equal(mcols(join_dist)$distance, mcols(join_test)$test)
  
  meta <- DataFrame('test' = 1:length(query))
  mcols(query) <- meta
  expect_error(join_fun(query, subject, distance = "test"), "columns already exist")
  
}

test_that("join_nearest behaves correctly on GRanges", {

  subject  <- data.frame(seqnames = "chr1",
               start = c(11,101),
               end = c(21, 200),
               name = c("a1", "a2"),
               strand = c("+", "-"),
               score = c(1,2)) %>%
           as_granges()
  query <- data.frame(seqnames = "chr1",
                        strand = c("+", "-", "+", "-"),
                        start = c(21,91,101,201),
                        end = c(30,101,110,210),
                        name = paste0("b", 1:4),
                        score = 1:4) %>%
                     as_granges()
  
  test_nearest_join_distance(query, subject, join_nearest)
  test_nearest_join_distance(query, subject, join_nearest_left)
  test_nearest_join_distance(query, subject, join_nearest_right)
  test_nearest_join_distance(query, subject, join_nearest_upstream)
  test_nearest_join_distance(query, subject, join_nearest_downstream)
  
})

test_that("join_nearest works on IRanges", {
  query_ir <- IRanges(c(6, 11, 1, 13, 18), width=c(2, 2, 2, 2, 2))
  subject_ir <- IRanges(c(1, 2, 9, 15, 15), width=c(4, 3, 2, 2, 3))
  
  test_nearest_join_distance(query_ir, subject_ir, join_nearest, ranges)
  test_nearest_join_distance(query_ir, subject_ir, join_nearest_left, ranges)
  test_nearest_join_distance(query_ir, subject_ir, join_nearest_right, ranges)
})

test_add_nearest_distance <- function(query, subject){
  dist <- add_nearest_distance(query, subject)
  test <- add_nearest_distance(query, subject, name = "test")
  
  meta <- DataFrame('test' = 1:length(query))
  mcols(query) <- meta
  
  # test errors on existing column
  expect_error(add_nearest_distance(query, subject, name = "test"), "already exists in destination")
  
  expect_equal(mcols(dist)$distance, mcols(test)$test)
}

test_that("add_nearest_distance works", {
  subject  <- data.frame(seqnames = "chr1",
               start = c(11,101),
               end = c(21, 200),
               name = c("a1", "a2"),
               strand = c("+", "-"),
               score = c(1,2)) %>%
           as_granges()
  query <- data.frame(seqnames = "chr1",
                        strand = c("+", "-", "+", "-"),
                        start = c(21,91,101,201),
                        end = c(30,101,110,210),
                        name = paste0("b", 1:4),
                        score = 1:4) %>%
                     as_granges()
  
  query_ir <- IRanges(c(6, 11, 1, 13, 18), width=c(2, 2, 2, 2, 2))
  subject_ir <- IRanges(c(1, 2, 9, 15, 15), width=c(4, 3, 2, 2, 3))

  test_add_nearest_distance(query, subject)
  test_add_nearest_distance(query_ir, subject_ir)
})
