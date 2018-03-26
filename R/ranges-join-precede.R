# ranges-join-precede.R

#' Find preceding Ranges
#'
#' @param x,y Ranges objects, which ranges in x precede those in y.
#' @param suffix A character vector of length two used to identify
#' metadata columns coming from x and y.
#'
#' @details By default `join_precede` will return the ranges
#' in x that come before the ranges in y and ignore any strand information.
#' The function `join_precede_right` will find all ranges in y
#' that are on the right-hand side of the ranges in x ignoring any strand
#' information. Finally, `join_precede_downstream` will find all ranges in y
#' that are that are downstream of the ranges in x. On the positive strand this
#' will result in ranges in y that are right of those in x and on the negative
#' strand it will result in ranges in y that are left of those in x.
#'
#' @return A Ranges object corresponding to the ranges in `y` that are
#' preceded by the ranges in `x`, all metadata is copied over from the
#' right-hand side ranges `y`.
#'
#' @examples
#' subject <- data.frame(start = c(5,10, 15,20), width = 5, gc = runif(4)) %>%
#'              as_iranges()
#' query <- data.frame(start = 2:6, width = 3:7, label = letters[1:5]) %>%
#'              as_iranges()
#'
#' join_precede(query, subject)
#'
#' query  <- data.frame(seqnames = "chr1",
#'                start = c(11,101),
#'                end = c(21, 200),
#'                name = c("a1", "a2"),
#'                strand = c("+", "-"),
#'                score = c(1,2)) %>%
#'            as_granges()
#' subject <- data.frame(seqnames = "chr1",
#'                       strand = c("+", "-", "+", "-"),
#'                       start = c(21,91,101,201),
#'                       end = c(30,101,110,210),
#'                       name = paste0("b", 1:4),
#'                       score = 1:4) %>%
#'                    as_granges()
#'
#' join_precede(query, subject)
#' join_precede_right(query, subject)
#' join_precede_downstream(query, subject)
#' @rdname ranges-precede
#' @importFrom IRanges precede
#' @export
join_precede <- function(x,y, suffix = c(".x", ".y")) { UseMethod("join_precede") }

#' @export
join_precede.IntegerRanges <- function(x,y, suffix = c(".x", ".y")) {
  hits <- precede(x,y)
  no_hits_id <- !is.na(hits)
  left <- x[no_hits_id, ]
  right <- y[hits[no_hits_id], ]
  mcols(left) <- mcols_overlaps_update(left, right, suffix)
  return(left)
}

#' @export
join_precede.GenomicRanges <- function(x,y, suffix = c(".x", ".y")) {
  hits <- precede(x,y, ignore.strand = TRUE)
  no_hits_id <- !is.na(hits)
  left <- x[no_hits_id, ]
  right <- y[hits[no_hits_id], ]
  mcols(left) <- mcols_overlaps_update(left, right, suffix)
  return(left)
}


#' @rdname ranges-precede
#' @importFrom IRanges precede
#' @export
join_precede_right <- function(x,y, suffix = c(".x", ".y")) { UseMethod("join_precede_right") }

#' @export
join_precede_right.IntegerRanges <- function(x,y, suffix = c(".x", ".y")) {
  hits <- precede(x,y, select = "all")
  left <- x[queryHits(hits), ]
  right <- y[subjectHits(hits), ]
  mcols(left) <- mcols_overlaps_update(left, right, suffix)
  left
}

#' @export
join_precede_right.GenomicRanges <- function(x,y, suffix = c(".x", ".y")) {
  hits <- precede(x,y, select = "all", ignore.strand = TRUE)
  left <- x[queryHits(hits), ]
  right <- y[subjectHits(hits), ]
  mcols(left) <- mcols_overlaps_update(left, right, suffix)
  left
}


#' @rdname ranges-precede
#' @importFrom IRanges precede
#' @export
join_precede_downstream <- function(x,y, suffix = c(".x", ".y")) {UseMethod("join_precede_downstream")}

#' @export
join_precede_downstream.GenomicRanges <- function(x,y, suffix = c(".x", ".y")) {
  hits <- precede(x,y, select = "all", ignore.strand = FALSE)
  left <- x[queryHits(hits), ]
  right <- y[subjectHits(hits), ]
  mcols(left) <- mcols_overlaps_update(left, right, suffix)
  left
}
