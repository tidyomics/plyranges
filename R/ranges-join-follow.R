# ranges-join-follow.R

#' Find following Ranges
#'
#' @param x,y Ranges objects, which ranges in x follow those in y.
#' @param suffix A character vector of length two used to identify
#' metadata columns coming from x and y.
#'
#' @details By default `join_follow` will find abritrary ranges
#' in y that are followed by ranges in x and ignore any strand information.
#' On the other hand `join_follow_left` will find all ranges in y
#' that are on the left-hand side of the ranges in x ignoring any strand
#' information. Finally, `join_follow_upstream` will find all ranges in x
#' that are that are upstream of the ranges in y. On the positive strand this
#' will result in ranges in y that are left of those in x and on the negative
#' strand it will result in ranges in y that are right of those in x.
#'
#' @return A Ranges object corresponding to the ranges in `x`` that are
#' followed by the ranges in `y`, all metadata is copied over from the
#' right-hand side ranges `y`.
#'
#' @rdname ranges-follow
#' @examples
#' query <- data.frame(start = c(5,10, 15,20), width = 5, gc = runif(4)) %>%
#'              as_iranges()
#' subject <- data.frame(start = 2:6, width = 3:7, label = letters[1:5]) %>%
#'              as_iranges()
#'
#' join_follow(query, subject)
#'
#' subject  <- data.frame(seqnames = "chr1",
#'                start = c(11,101),
#'                end = c(21, 200),
#'                name = c("a1", "a2"),
#'                strand = c("+", "-"),
#'                score = c(1,2)) %>%
#'            as_granges()
#' query <- data.frame(seqnames = "chr1",
#'                       strand = c("+", "-", "+", "-"),
#'                       start = c(21,91,101,201),
#'                       end = c(30,101,110,210),
#'                       name = paste0("b", 1:4),
#'                       score = 1:4) %>%
#'                    as_granges()
#'
#' join_follow(query, subject)
#' join_follow_left(query, subject)
#' join_follow_upstream(query, subject)
#'
#' @importFrom IRanges follow
#' @export
join_follow <- function(x,y, suffix = c(".x", ".y")) { UseMethod("join_follow") }

#' @export
join_follow.IntegerRanges <- function(x,y, suffix = c(".x", ".y")) {
  hits <- follow(x,y)
  expand_by_hits(x, y, suffix, hits)
}

#' @export
join_follow.GenomicRanges <- function(x,y, suffix = c(".x", ".y")) {
  hits <- follow(x,y, ignore.strand = TRUE)
  expand_by_hits(x, y, suffix, hits)
}

#' @rdname ranges-follow
#' @importFrom IRanges follow
#' @export
join_follow_left <- function(x,y, suffix = c(".x", ".y")) { UseMethod("join_follow_left") }

#' @export
join_follow_left.IntegerRanges <- function(x,y, suffix = c(".x", ".y")) {
  hits <- follow(x,y, select = "all")
  expand_by_hits(x,y, suffix, hits)
}


#' @export
join_follow_left.GenomicRanges <- function(x,y, suffix = c(".x", ".y")) {
  hits <- follow(x,y, select = "all", ignore.strand = TRUE)
  expand_by_hits(x,y,suffix, hits)
}

#' @rdname ranges-follow
#' @importFrom IRanges follow
#' @export
join_follow_upstream <- function(x,y, suffix = c(".x", ".y")) {UseMethod("join_follow_upstream")}

#' @export
join_follow_upstream.GenomicRanges <- function(x,y, suffix = c(".x", ".y")) {
  hits <- follow(x,y, select = "all", ignore.strand = FALSE)
  expand_by_hits(x,y,suffix, hits)
}
