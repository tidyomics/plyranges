# ranges-anchors-class.R

# A virtual class that delegates anchoring, so we can decorate
# a Ranges derivative object with an anchor.
validAnchoredGenomicRanges <- function(object) {
  valid_anchors <- c("start", "end", "center", "centre", "5p", "3p")
  
  if (length(object@anchor) > 1L) {
    paste("anchor must be character vector of length 1.")
  }
  
  if (!(object@anchor %in% valid_anchors)) {
    paste(object@anchor, "is not a valid anchor.")
  }
}

setClass("AnchoredGenomicRanges",
         representation = representation(anchor = "character"),
         contains = c("VIRTUAL", "DelegatingGenomicRanges"),
         validity = validAnchoredGenomicRanges)


setMethod("seqnames", "AnchoredGenomicRanges",
          function(x) seqnames(x@delegate))

setMethod("ranges", "AnchoredGenomicRanges",
          function(x, ...) ranges(x@delegate, ...))

setMethod("strand", "AnchoredGenomicRanges", function(x) strand(x@delegate))
setMethod("seqinfo", "AnchoredGenomicRanges", function(x) seqinfo(x@delegate))
setMethod("mcols", "AnchoredGenomicRanges", function(x, ...) {
  mcols(x@delegate, ...)
})
setMethod("update", "AnchoredGenomicRanges", function (object, ...) {
  object@delegate <- update(object@delegate, ...)
  object
})

setMethod("show", "AnchoredGenomicRanges", function(object) {
  output <- c("", utils::capture.output(show(object@delegate)))
  output[1] <- output[2]
  anchor <- object@anchor
  output[2] <- paste("Anchored by:", anchor)
  cat(output, sep = "\n")
})

new_anchored_gr <- function(rng, anchor) {
  cls <- class(rng)
  new_cls <- paste0(cls, "Anchored")
  setClass(new_cls, contains =  "AnchoredGenomicRanges")
  new(new_cls, delegate = rng, anchor = anchor)
}

# equivalent of DelegatingGenomicRanges for IRanges objects
setClass("DelegatingIntegerRanges",
         representation = representation(delegate = "IntegerRanges"),
         contains = c("VIRTUAL", "IntegerRanges"))

validAnchoredIntegerRanges <- function(object) {
  valid_anchors <- c("start", "end", "center", "centre")
  
  if (length(object@anchor) > 1L) {
    paste("anchor must be character vector of length 1.")
  }
  
  if (!(object@anchor %in% valid_anchors)) {
    paste(object@anchor, "is not a valid anchor.")
  }
}

setClass("AnchoredIntegerRanges",
         slot = c(anchor = "character"),
         contains = c("VIRTUAL", "IntegerRanges"),
         validity = validAnchoredIntegerRanges)

setMethod("start", "AnchoredIntegerRanges", function(x, ...) start(x@delegate))
setMethod("end", "AnchoredIntegerRanges", function(x, ...) end(x@delegate))
setMethod("width", "AnchoredIntegerRange", function(x) width(x@delegate))

