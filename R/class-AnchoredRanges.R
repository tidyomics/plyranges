# A virtual class that delegates anchoring, so we can decorate
# a Ranges derivative object with an anchor.
# --- GenomicRanges and friends ---
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
         contains = "DelegatingGenomicRanges",
         validity = validAnchoredGenomicRanges)

# constructor
new_anchored_gr <- function(rng, anchor) {
  new("AnchoredGenomicRanges", 
      elementMetadata =  S4Vectors:::make_zero_col_DataFrame(length(rng)), 
      delegate = rng, 
      anchor = anchor)
}

# mcols method for DelegatingGenomicRanges
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

# --- IntegerRanges and friends ---
# Equivalent of DelegatingGenomicRanges for IntegerRanges class
setClass("DelegatingIntegerRanges",
          representation = representation(delegate = "IntegerRanges"),
          contains = c("VIRTUAL", "IntegerRanges"))

# methods for DelegatingIntegerRanges
setMethod("start", "DelegatingIntegerRanges", function(x, ...) start(x@delegate))
setMethod("end", "DelegatingIntegerRanges", function(x, ...) end(x@delegate))
setMethod("width", "DelegatingIntegerRanges", function(x) width(x@delegate))
setMethod("mcols", "DelegatingIntegerRanges", function(x, ...) {
    mcols(x@delegate, ...)
})

# --- ready for our AnchoredIntegerRanges ---
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



setMethod("show", "AnchoredIntegerRanges", function(object) {
  output <- c("", utils::capture.output(show(object@delegate)))
  output[1] <- output[2]
  anchor <- object@anchor
  output[2] <- paste("Anchored by:", anchor)
  cat(output, sep = "\n")
})

new_anchored_ir <- function(rng, anchor) {
  new("AnchoredIntegerRanges", delegate = rng, anchor = anchor)
}