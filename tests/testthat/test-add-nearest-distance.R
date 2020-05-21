compare_dist_vals <- function(add, join){
  # Checks that joins and adds produce the same values (minus NA's)
  add_dist <- mcols(join)$distance
  add_dist <- add_dist[!is.na(add_dist)]
  join_dist <- mcols(join)$distance
  
  expect_true(identical(add_dist, join_dist))
}

test_that("add_nearest_distance works on GRanges", {
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
  
  near <- add_nearest_distance(query, subject)
  near_left <- add_nearest_distance_left(query, subject)
  near_right <- add_nearest_distance_right(query, subject)
  near_up <- add_nearest_distance_upstream(query, subject)
  near_down <- add_nearest_distance_downstream(query, subject)
  
  j_near <- join_nearest(query, subject, distance = TRUE)
  j_near_left <- join_nearest_left(query, subject, distance = TRUE)
  j_near_right <- join_nearest_right(query, subject, distance = TRUE)
  j_near_up <- join_nearest_upstream(query, subject, distance = TRUE)
  j_near_down <- join_nearest_downstream(query, subject, distance = TRUE)
  
  
  compare_dist_vals(near, j_near)
  compare_dist_vals(near_left, j_near_left)
  compare_dist_vals(near_right, j_near_right)
  compare_dist_vals(near_up, j_near_up)
  compare_dist_vals(near_down, j_near_down)
  
})

test_that("add_nearest_distance works on IRanges", {
  query_ir <- IRanges(c(6, 11, 1, 13, 18), width=c(2, 2, 2, 2, 2))
  subject_ir <- IRanges(c(1, 2, 9, 15, 15), width=c(4, 3, 2, 2, 3))
  
  near <- add_nearest_distance(query_ir, subject_ir)
  near_left <- add_nearest_distance_left(query_ir, subject_ir)
  near_right <- add_nearest_distance_right(query_ir, subject_ir)
  
  j_near <- join_nearest(query_ir, subject_ir, distance = TRUE)
  j_near_left <- join_nearest_left(query_ir, subject_ir, distance = TRUE)
  j_near_right <- join_nearest_right(query_ir, subject_ir, distance = TRUE)
  
  compare_dist_vals(near, j_near)
  compare_dist_vals(near_left, j_near_left)
  compare_dist_vals(near_right, j_near_right)
 
  # no up/down method for IRanges 
  expect_error(add_nearest_distance_upstream(query_ir, subject_ir), "no applicable method")
  expect_error(add_nearest_distance_downstream(query_ir, subject_ir), "no applicable method")
  
})
