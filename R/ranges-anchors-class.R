# ranges-anchors-class.R

# A virtual class that delegates anchoring, so we can decorate
# a Ranges derivative object with an anchor.
# -- for GenomicRanges derivatives
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
         contains = c("DelegatingGenomicRanges"),
         validity = validAnchoredGenomicRanges)

setMethod("seqnames", "AnchoredGenomicRanges",
          function(x) seqnames(x@delegate))

setMethod("ranges", "AnchoredGenomicRanges",
          function(x, ...) ranges(x@delegate, ...))

setMethod("strand", "AnchoredGenomicRanges", function(x) strand(x@delegate))

setMethod("seqinfo", "AnchoredGenomicRanges", function(x) seqinfo(x@delegate))

setMethod("mcols", "DelegatingGenomicRanges", function(x, ...) {
  mcols(x@delegate, ...)
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
  setClass(new_cls, contains = c("AnchoredGenomicRanges"), where = parent.frame())
  new(new_cls, 
      elementMetadata =  S4Vectors:::make_zero_col_DataFrame(length(rng)), 
      delegate = rng, 
      anchor = anchor)
}

# -- for IntegerRanges derivatives
# # equivalent of DelegatingGenomicRanges for IRanges objects
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
         contains = "DelegatingIntegerRanges",
         validity = validAnchoredIntegerRanges)

setMethod("start", "AnchoredIntegerRanges", function(x, ...) start(x@delegate))
setMethod("end", "AnchoredIntegerRanges", function(x, ...) end(x@delegate))
setMethod("width", "AnchoredIntegerRanges", function(x) width(x@delegate))

setMethod("show", "AnchoredIntegerRanges", function(object) {
  output <- c("", utils::capture.output(show(object@delegate)))
  output[1] <- output[2]
  anchor <- object@anchor
  output[2] <- paste("Anchored by:", anchor)
  cat(output, sep = "\n")
})

new_anchored_ir <- function(rng, anchor) {
  cls <- class(rng)
  new_cls <- paste0(cls, "Anchored")
  setClass(new_cls, contains =  "AnchoredIntegerRanges", where = parent.frame())
  new(new_cls, delegate = rng, anchor = anchor)
}