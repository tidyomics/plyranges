context("reduce_ranges method")

test_that("matches IRanges/GenomicRanges", {
  x <- IRanges()
  expect_identical(x, reduce_ranges(x))

  x <- IRanges(c(1:4, 10:11, 11), width=c(0,1,1,0,0,0,1))
  mcols(x) <- DataFrame(mapping = paste0("a", seq_along(x)))
  target <- IRanges(c(1:2, 10:11), width=c(0,2,0,1))
  mcols(target) <- DataFrame(mapping=CharacterList("a1", c("a2","a3","a4"), "a5", c("a6","a7")))
  expect_identical(reduce_ranges(x, mapping = mapping),
                   target)

  # drop.empty.ranges is just a filter
  current <- x %>%
    filter(width > 0) %>%
    reduce_ranges(mapping = mapping)
  target <- reduce(x, drop.empty.ranges=TRUE)
  mcols(target) <- DataFrame(mapping = CharacterList(c("a2","a3"), "a7"))

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

  mcols(target)$mapping <- CharacterList("a6", "a1", "a5", c("a2","a3"), 
                                         "a4", c("a7","a8"), c("a9","a10"))

  expect_identical(reduce_ranges_directed(gr, mapping = name),
                   target)
  target <- GRanges(Rle(c("chr1", "chr2", "chr3"), c(1,1,1)),
                    IRanges(start = c(1,2,7), end = c(10,10,10)))
  expect_identical(reduce_ranges(gr), target)

  mcols(target)$mapping <- CharacterList(
    c("a1","a5","a6"), 
    c("a2","a3","a4"), 
    c("a7","a8","a9","a10"))

  expect_identical(reduce_ranges(gr, mapping = name),
                   target)


})


test_that("non-standard evaluation works as expected",{
  oldwd <- getwd()
  setwd(system.file("unitTests", "data", "merge", package="HelloRanges"))

  a <- read_bed("a.bed")

  exp <- reduce(a, ignore.strand = TRUE)

  expect_identical(exp, reduce_ranges(a))

  exp <- reduce(a, with.revmap=TRUE, ignore.strand = TRUE)
  mcols(exp)[["seqnames.count"]] <- lengths(exp$revmap)
  exp <- exp %>% select(-revmap)

  expect_identical(exp, a %>% reduce_ranges(seqnames.count = n()))

  mcols(a)$name <- paste0("a", 1:4)
  exp <- reduce(a, with.revmap=TRUE)
  mcols(exp)[["name.collapse"]] <- CharacterList("a1", c("a2","a3","a4"))
  exp <- exp %>% select(-revmap)

  expect_identical(exp, a %>%
                     reduce_ranges(name.collapse = name)
                   )
  a <- read_bed("a.full.bed")
  exp <- reduce(a, with.revmap=TRUE, ignore.strand=TRUE)
  mcols(exp)[["name.collapse"]] <- CharacterList("a1", c("a2", "a3","a4"), "a1", "a2", c("a3","a4"))
  mcols(exp)[["score.sum"]] <- c(1,9,5,6,15)

  exp <- exp %>% select(-revmap)

  expect_identical(exp, reduce_ranges(a,
                                      name.collapse = name,
                                      score.sum = sum(score)))

  exp <- reduce(a, with.revmap=TRUE, ignore.strand=TRUE)
  mcols(exp)[["score.count"]] <- c(1L,3L,1L,1L,2L)
  mcols(exp)[["score.sum"]] <- c(1,9,5,6,15)
  exp <- exp %>% select(-revmap)

  expect_identical(exp, reduce_ranges(a, score.count = n(),
                                      score.sum = sum(score)))
  setwd(oldwd)
})

test_that("grouping then reducing works as expected", {
  oldwd <- getwd()
  setwd(system.file("unitTests", "data", "multiinter", package="HelloRanges"))
  bed_files <- list.files(pattern = ".bed$")
  # GRangesList
  gr_l <- as(lapply(bed_files, function(x) {
    mutate(read_bed(x), grp = sub(".bed$", "", basename(x)))
  }), "GRangesList")
  names(gr_l) <- sub(".bed$", "", basename(bed_files))
  gr_l_reduced <- reduce(gr_l)
  correct_gr <- IRanges::stack(gr_l_reduced, "grp") %>%
    mutate(grp = as.character(grp))

  # GroupedGRanges
  gr_by_group <- unlist(gr_l, use.names = FALSE) %>% group_by(grp)
  test_gr <- reduce_ranges(gr_by_group)

  expect_identical(correct_gr, test_gr)

  # with an operation matches  revmap length
  gr_l_reduced <- reduce(gr_l, with.revmap = TRUE)
  correct_n <- lengths(IRanges::stack(gr_l_reduced)$revmap)
  test_n <- reduce_ranges(gr_by_group, n = n())$n
  expect_identical(correct_n, test_n)
  setwd(oldwd)
  
})


test_that("expected behaviour for grouped filter w reduce #37",
          {
            # see https://github.com/sa-lee/plyranges/issues/37
            set.seed(2019)
            n <- 10
            r <- GRanges(seqnames = rep("chr1", n),
                         ranges = IRanges(start = sample(20, n, replace = TRUE),
                                          width = sample(6,  n, replace = TRUE))
            )
            mcols(r) <- data.frame(score = runif(n, 0, 100), 
                                   condition = rep_len(c("One","Two"), n))
            red1 <- r %>% group_by(condition) %>% reduce_ranges()
            
            exp <- S4Vectors::split(r, r$condition) %>% 
              reduce() %>%
              unlist()
            exp <- exp %>%
              mutate(condition = as.factor(names(.)))
            names(exp) <- NULL
            
            expect_identical(red1, exp)
            
            red2 <- r %>%  
              group_by(condition) %>% 
              filter(score > 2) %>% 
              reduce_ranges()
            
            exp <- S4Vectors::split(r, r$condition) 
            is_gt2 <- as(lapply(exp, function(x) x$score > 2), "List")
            exp <- exp[is_gt2] %>% reduce() %>% unlist()
            exp <- exp %>%
              mutate(condition = as.factor(names(.)))
            names(exp) <- NULL
            
            expect_identical(red2, exp)
          }
)