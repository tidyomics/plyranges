# A virtual class that delegates anchoring, so we can decorate
# a Ranges derivative object with an anchor.
# --- GenomicRanges and friends ---
validAnchoredRanges <- function(object, valid_anchors =  c("start", "end", "center", "centre", "5p", "3p")) {  
  if (length(object@anchor) > 1L) {
    paste("anchor must be character vector of length 1.")
  }
  
  if (!(object@anchor %in% valid_anchors)) {
    paste(object@anchor, "is not a valid anchor.")
  }
}

setClass("AnchoredGenomicRanges",
         representation = representation(anchor = "character"),
         contains = "DelegatingGenomicRanges")

setValidity("AnchoredGenomicRanges", function(object) {
  validAnchoredRanges(object)
})

# constructor
initialize_AnchoredRanges <- function(.Object, delegate, anchor) {
  .Object@delegate <- delegate
  .Object@anchor <- anchor
  .Object
}

setMethod("initialize", "AnchoredGenomicRanges", 
          function(.Object, delegate, anchor, ...) {
            initialize_AnchoredRanges(.Object, delegate, anchor)
})

new_anchored_gr <- function(rng, anchor) {
  new("AnchoredGenomicRanges", delegate = rng, anchor = anchor)
}

# mcols method for DelegatingGenomicRanges
setMethod("mcols", "DelegatingGenomicRanges", function(x, ...) {
  mcols(x@delegate, ...)
})

show_AnchoredRanges <- function(object) {
  output <- c("", utils::capture.output(show(object@delegate)))
  output[1] <- output[2]
  anchor <- object@anchor
  output[2] <- paste("Anchored by:", anchor)
  cat(output, sep = "\n")
}
setMethod("show", "AnchoredGenomicRanges", show_AnchoredRanges)

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

#' @importFrom methods callNextMethod
setMethod("mcols<-", "DelegatingGenomicRanges", function(x, ..., value) {
  x@delegate <- callNextMethod(x = x@delegate, ..., value = value)
  x
})


setMethod("mcols<-", "DelegatingIntegerRanges", function(x, ..., value) {
  x@delegate <- callNextMethod(x = x@delegate, ..., value = value)
  x
})

#' @importFrom S4Vectors parallelVectorNames
setMethod("parallelVectorNames", "DelegatingGenomicRanges", function(x) {
  callNextMethod(x = x@delegate)
})

#' @importFrom S4Vectors parallelVectorNames
setMethod("parallelVectorNames", "DelegatingIntegerRanges", function(x) {
  callNextMethod(x = x@delegate)
})

# --- ready for our AnchoredIntegerRanges ---
setClass("AnchoredIntegerRanges",
         slot = c(anchor = "character"),
         contains = "DelegatingIntegerRanges")

setValidity("AnchoredIntegerRanges", function(object) {
  validAnchoredRanges(object, c("start", "end", "center", "centre"))
})

setMethod("initialize", "AnchoredIntegerRanges", 
          function(.Object, delegate, anchor, ...) {
            initialize_AnchoredRanges(.Object, delegate, anchor)
})

setMethod("show", "AnchoredIntegerRanges", show_AnchoredRanges)

new_anchored_ir <- function(rng, anchor) {
  new("AnchoredIntegerRanges", delegate = rng, anchor = anchor)
}