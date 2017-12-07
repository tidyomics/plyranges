# ranges-anchors.R

#' An S4 class to represent anchored GRanges
#'
#' @slot anchor a character(1) containing anchor. \code{anchor} must be one of
#' \code{c("start", "end", "center", "3p", "5p")}.
#' @importFrom methods setClass setValidity setMethod show
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

#' An S4 class to represent anchored IRanges
#'
#' @slot anchor a character(1) containing anchor. \code{anchor} must be one of
#' \code{c("start", "end", "center")}.
#' @rdname IRangesAnchored-class
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

#' Fixing a Range by coordinates or strand
#' @param x a Ranges object
#'
#' @details Anchoring will fix a Ranges start, end, or center positions,
#' so these positions will remain the same when performing arithimetic.
#' Similiarly, (\code{anchor_3p}) will fix the start
#' for the negative strand, while \code{anchor_5p} will fix the end for the
#' positive strand. To see what has been anchored use
#' \code{anchor}.
#' @return A RangesAnchored object
#' @importFrom methods is new
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
