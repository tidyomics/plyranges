# ranges-anchors.R
#' @importFrom methods is
is_ranges <- function(x) {
  is(x, "Ranges") || is(x, "GenomicRanges")
}


#' Fixing a Range by coordinates or strand
#' @param x a Ranges object
#' @details Anchoring will fix a Ranges start, end, or center positions,
#' so these positions will remain the same when performing arithimetic.
#' Similiarly, fixing by positive or negative strand will leave the intervals
#' with those strand values unchanged. To see what has been anchored use
#' \code{anchors}.
#' @rdname anchor
#' @export
anchors <- function(x) {
  stopifnot(is_ranges(x))
  metadata(x)[["anchor"]]
}

#' @describeIn anchor Fix by start position
#' @export
anchor_start <- function(x) {
  stopifnot(is_ranges(x))
  metadata(x)$anchor <- append(metadata(x)$anchor, "start")
  x
}

#' @describeIn anchor Fix by end position
#' @export
anchor_end <- function(x) {
  stopifnot(is_ranges(x))
  metadata(x)$anchor <- append(metadata(x)$anchor, "end")
  x
}

#' @describeIn anchor Fix by center position
#' @export
anchor_center <- function(x) {
  stopifnot(is_ranges(x))
  metadata(x)$anchor <- append(metadata(x)$anchor, "center")
  x
}

#' @describeIn anchor Fix by positive strand
#' @export
anchor_3p <- function(x) {
  stopifnot(is_ranges(x))
  metadata(x)$anchor <- append(metadata(x)$anchor, "3p")
  x
}

#' @describeIn anchor Fix by negative strand
#' @export
anchor_5p <- function(x) {
  stopifnot(is_ranges(x))
  metadata(x)$anchor <- append(metadata(x)$anchor, "5p")
  x
}
