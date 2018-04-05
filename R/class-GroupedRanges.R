validGroupedGenomicRanges <- function(object) {
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

#' @rdname group_by-ranges
#' @importFrom IRanges extractList splitAsList
#' @importFrom BiocGenerics unlist
#' @export
setClass("GroupedGenomicRanges",
         slot = c(groups = "list",
                  inx = "IntegerList"),
         contains = "DelegatingGenomicRanges",
         validity = validGroupedGenomicRanges)

# generates index for grouping variables
make_group_inx <- function(rng, ...) {
    capture_groups <- quos(...)
    group_names <- unlist(lapply(capture_groups, quo_name))
    groups <- rlang::syms(group_names)
    names(groups) <- group_names
    # eval groups
    os <- overscope_ranges(rng)
    on.exit(rlang::overscope_clean(os))
    groups_values <- as(overscope_eval_update(os, groups), "DataFrame")
    inx <- IRanges::splitAsList(seq_along(rng), groups_values, drop = TRUE)
    mcols(inx) <- BiocGenerics::unlist(extractList(groups_values, 
                              endoapply(inx, function(.) .[1])))
    return(list(groups = groups, inx = inx))
}

# constructor (passed to `group_by`)
new_grouped_gr <- function(rng, ...) {
    groupings <- make_group_inx(rng, ...)
    # instantiate class
    new("GroupedGenomicRanges",
        elementMetadata =  S4Vectors:::make_zero_col_DataFrame(length(rng)),
        delegate = rng, 
        groups = groupings$groups,
        inx = groupings$inx)
}

show_GroupedRanges <- function(object) {
    groups <- unlist(lapply(object@groups, as.character))
    groups <- paste(groups, collapse = ", ")
    output <- c("", utils::capture.output(show(object@delegate)))
    output[1] <- output[2]
    output[2] <- paste("Groups:", groups, paste0("[", length(object@inx), "]"))
    cat(output, sep = "\n")
}
setMethod("show", "GroupedGenomicRanges", function(object) { 
    show_GroupedRanges(object)
})

# --- GroupedIntegerRanges ---
validGroupedIntegerRanges <- function(object) {
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

#' @rdname group_by-ranges
#' @export
setClass("GroupedIntegerRanges",
         slot = c(groups = "list",
                  inx = "IntegerList"),
         contains = "DelegatingIntegerRanges",
         validity = validGroupedIntegerRanges)

new_grouped_ir <- function(rng, ...) {
    groupings <- make_group_inx(rng, ...)
    # instantiate class
    new("GroupedIntegerRanges",
        delegate = rng, 
        groups = groupings$groups,
        inx = groupings$inx)
}

setMethod("show", "GroupedIntegerRanges", function(object) {
    show_GroupedRanges(object)
})
