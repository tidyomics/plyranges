# ranges-anchors.R

#' Anchored Ranges objects
#'
#' @description The \code{GRangesAnchored} class and the \code{IRangesAnchored}
#' class allow components of a \code{GRanges} or \code{IRanges} (start, end, center)
#' to be held fixed.
#'
#' @details Anchoring will fix a Ranges start, end, or center positions,
#' so these positions will remain the same when performing arithimetic.
#' For \code{GRanges} objects, the function
#' (\code{anchor_3p}) will fix the start for the negative strand,
#' while \code{anchor_5p} will fix the end for the
#' positive strand. Anchoring modifies how arithmetic is performed, for example
#' modifying the width of a range with \code{set_width} or stretching a
#' range with \code{stretch}.
#'
#' @param x a Ranges object
#'
#' @section Constructors:
#' Depending on how you want to fix the components of a Ranges, there are
#' five ways to construct a RangesAnchored class. Here \code{x} is either
#' an \code{IRanges} or \code{GRanges} object.
#' \itemize{
#'    \item{\code{anchor_start(x)}}{Fix the start coordinates}
#'    \item{\code{anchor_end(x)}}{Fix the end coordinates}
#'    \item{\code{anchor_center(x)}}{Fix the center coordinates}
#'    \item{\code{anchor_3p(x)}}{On the negative strand fix the start coordinates,
#'    and for positive or unstranded ranges fix the end coordinates.}
#'    \item{\code{anchor_5p(x)}}{On the positive or unstranded ranges fix the start coordinates,
#'    coordinates and for negative stranded ranges fix the end coordinates.}
#' }
#'
#' @section Accessors:
#' To see what has been anchored use the function \code{anchor}.
#' This will return a character vector containing a valid anchor.
#' It will be set to one of \code{c("start", "end", "center")} for an
#' \code{IRanges} object or one of
#' \code{c("start", "end", "center", "3p", "5p")} for a \code{GRanges} object.
#'
#' @seealso \link{set_width}, \link{stretch}
#'
#' @examples
#' df <- data.frame(start = 1:10, width = 5)
#' rng <- as_iranges(df)
#' rng_by_start <- anchor_start(rng)
#' rng_by_start
#' anchor(rng_by_start)
#' set_width(rng_by_start, 3L)
#' grng <- as_granges(df,
#'                    seqnames = "chr1",
#'                    strand = c(rep("-", 5), rep("+", 5)))
#' rng_by_5p <- anchor_5p(grng)
#' rng_by_5p
#' set_width(rng_by_5p, 3L)
#'
#' @importFrom methods setClass setValidity setMethod show
#' @rdname ranges-anchor
#' @export
setClass("GRangesAnchored",
         slot = c(anchor = "character"),
         contains = "GRanges")

validGRangesAnchored <- function(object) {
  valid_anchors <- c("start", "end", "center", "centre", "5p", "3p")

  if (length(object@anchor) > 1L) {
    paste("anchor must be character vector of length 1.")
  }

  if (!(object@anchor %in% valid_anchors)) {
    paste(object@anchor, "is not a valid anchor.")
  }
}


setValidity("GRangesAnchored", validGRangesAnchored)

setMethod("show", "GRangesAnchored", function(object) {
  output <- c("", utils::capture.output(show(as(object, "GenomicRanges"))))
  output[1] <- output[2]
  anchor <- object@anchor
  output[2] <- paste("Anchored by:", anchor)
  cat(output, sep = "\n")
})


#' @rdname ranges-anchor
#' @importFrom methods setClass setValidity setMethod
#' @export
setClass("IRangesAnchored",
         slot = c(anchor = "character"),
         contains = "IRanges")

validIRangesAnchored <- function(object) {
  valid_anchors <- c("start", "end", "center", "centre")

  if (length(object@anchor) > 1L) {
    paste("anchor must be character vector of length 1.")
  }

  if (!(object@anchor %in% valid_anchors)) {
    paste(object@anchor, "is not a valid anchor.")
  }
}

setValidity("IRangesAnchored", validIRangesAnchored)

setMethod("show", "IRangesAnchored", function(object) {
  output <- c("", utils::capture.output(show(as(object, "Ranges"))))
  output[1] <- output[2]
  anchor <- object@anchor
  output[2] <- paste("Anchored by:", anchor)
  cat(output, sep = "\n")
})


#' @rdname ranges-anchor
#' @export
anchor <- function(x) {
  if (!(is(x, "GRangesAnchored") | is(x, "IRangesAnchored"))) {
    return(NULL)
  } else {
    x@anchor
  }
}

#' @rdname ranges-anchor
#' @export
anchor_start <- function(x) { UseMethod("anchor_start") }

#' @export
anchor_start.Ranges <- function(x) {
  new("IRangesAnchored", anchor = "start", x)
}

#' @export
anchor_start.GenomicRanges <- function(x) {
  new("GRangesAnchored", anchor = "start", x)
}


#' @rdname ranges-anchor
#' @export
anchor_end <- function(x) { UseMethod("anchor_end") }

#' @export
anchor_end.Ranges <- function(x) {
  new("IRangesAnchored", anchor = "end", x)
}

#' @export
anchor_end.GenomicRanges <- function(x) {
  new("GRangesAnchored", anchor = "end", x)
}

#' @rdname ranges-anchor
#' @export
anchor_center <- function(x) { UseMethod("anchor_center") }

#' @export
anchor_center.Ranges <- function(x) {
  new("IRangesAnchored", anchor = "center", x)
}

#' @export
anchor_center.GenomicRanges <- function(x) {
  new("GRangesAnchored", anchor = "center", x)
}

#' @rdname ranges-anchor
#' @export
anchor_centre <- function(x) { UseMethod("anchor_center") }

#' @rdname ranges-anchor
#' @export
anchor_3p <- function(x) { UseMethod("anchor_3p") }

#' @export
anchor_3p.GenomicRanges <- function(x) {
  new("GRangesAnchored", anchor = "3p", x)
}

#' @rdname ranges-anchor
#' @export
anchor_5p <- function(x) { UseMethod("anchor_5p") }

#' @export
anchor_5p.GenomicRanges <- function(x) {
  new("GRangesAnchored", anchor = "5p", x)
}
