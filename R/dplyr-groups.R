

setClass("GroupedGRanges",
         slot = c(groups = "list"),
         contains = "GRanges")


validGroupedGRanges <- function(object) {
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

  TRUE

}


setValidity("GroupedGRanges", validGroupedGRanges)

groups.GroupedGRanges <- function(x) { x@groups }

setMethod("show", "GroupedGRanges", function(object) {
  groups <- unlist(lapply(groups(object), as.character))
  groups <- paste(groups, collapse = ", ")
  ranges_print <- c("", capture.output(GenomicRanges:::show_GenomicRanges(object,
                                 margin = "  ",
                                 print.classinfo = TRUE,
                                 print.seqinfo = TRUE)))
  ranges_print[1] <- ranges_print[2]
  ranges_print[2] <- paste("Groups:", groups)
  cat(ranges_print, sep = "\n")

})


#' @importFrom rlang quo_name quos syms
group_by.GRanges <- function(.data, ...) {
  capture_groups <- quos(...)
  groups <- lapply(capture_groups, function(x) quo_name(x))
  groups <- syms(groups)

  new("GroupedGRanges", .data, groups = groups)

}
