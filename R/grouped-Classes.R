
#' An S4 class to represent grouped GRanges
#'
#' @slot groups a list of names
#' @seealso \code{\link{group_by}}
#' @rdname GRangesGrouped-class
#' @importFrom methods setClass setValidity setMethod
#' @export
setClass("GRangesGrouped",
         slot = c(groups = "list"),
         contains = "GRanges")


validGRangesGrouped <- function(object) {
  check_valid_names <- all(unlist(lapply(object@groups, is.name)))

  if (!check_valid_names) {
    paste("Invalid groups slot: groups must be names")
  }

  group_names <- unlist(lapply(object@groups, as.character))
  check_valid_groups <- !(group_names %in% c(names(mcols(object)),
                                             "seqnames", "strand", "ranges",
                                             "start", "end", "width"))

  if (any(check_valid_groups)) {
    paste("Invalid groups slot:",
          paste(group_names[check_valid_groups], collapse = ","),
          "not found in data.")
  }

}


setValidity("GRangesGrouped", validGRangesGrouped)

setMethod("show", "GRangesGrouped", function(object) {
  groups <- unlist(lapply(object@groups, as.character))
  groups <- paste(groups, collapse = ", ")
  ranges_print <- c("", utils::capture.output(GenomicRanges:::show_GenomicRanges(object,
                                                                          margin = "  ",
                                                                          print.classinfo = TRUE,
                                                                          print.seqinfo = TRUE)))
  ranges_print[1] <- ranges_print[2]
  ranges_print[2] <- paste("Groups:", groups)
  cat(ranges_print, sep = "\n")

})

#' An S4 class to represent grouped IRanges
#'
#' @slot groups a list of names
#' @seealso \code{\link{group_by}}
#' @rdname IRangesGrouped-class
#' @export
setClass("IRangesGrouped",
         slot = c(groups = "list"),
         contains = "IRanges")


validIRangesGrouped <- function(object) {
  check_valid_names <- all(unlist(lapply(object@groups, is.name)))

  if (!check_valid_names) {
    paste("Invalid groups slot: groups must be names")
  }

  group_names <- unlist(lapply(object@groups, as.character))
  check_valid_groups <- !(group_names %in% c(names(mcols(object)),
                                             "start", "end", "width"))

  if (any(check_valid_groups)) {
    paste("Invalid groups slot:",
          paste(group_names[check_valid_groups], collapse = ","),
          "not found in data.")
  }
}

setValidity("IRangesGrouped", validIRangesGrouped)

setMethod("show", "IRangesGrouped", function(object) {
  groups <- unlist(lapply(object@groups, as.character))
  groups <- paste(groups, collapse = ", ")
  ranges_print <- c("", utils::capture.output(IRanges:::showRanges(object, margin = "  ",
                                                            print.classinfo = TRUE
  )))
  ranges_print[1] <- ranges_print[2]
  ranges_print[2] <- paste("Groups:", groups)
  cat(ranges_print, sep = "\n")

})

