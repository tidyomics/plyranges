validGroupedRanges <- function(object) {
  check_valid_names <- all(unlist(lapply(object@groups, is.name)))

  if (!check_valid_names) {
    paste("Invalid groups slot: groups must be names")
  }

  group_names <- unlist(lapply(object@groups, as.character))
  check_valid_groups <- !(group_names %in% ranges_vars(object))

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
         validity = validGroupedRanges)

initialize_GroupedRanges <- function(.Object, delegate, groups, inx) {
    stopifnot(is(delegate, "Ranges"))
    .Object@delegate <- delegate
    .Object@groups <- groups
    .Object@inx <- inx
    .Object
}

#' @importFrom methods setMethod initialize
setMethod("initialize", "GroupedGenomicRanges",
    function(.Object, delegate, groups, inx, ...) {
        initialize_GroupedRanges(.Object, delegate, groups, inx)
    }
)

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

# --- group-by backend ---
# generates index for grouping variables
make_group_inx <- function(rng, ...) {
    capture_groups <- rlang::enquos(...)
    group_names <- unlist(lapply(capture_groups, quo_name))
    groups <- rlang::syms(group_names)
    capture_groups <- lapply(groups, rlang::as_quosure)
    names(capture_groups) <- group_names
    # eval groups
    os <- overscope_ranges(rng)
    groups_values <- as(lapply(capture_groups, eval_tidy, data = os), 
                        "DataFrame")
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
        delegate = rng, 
        groups = groupings$groups,
        inx = groupings$inx)
}


# --- GroupedIntegerRanges ---

#' @rdname group_by-ranges
#' @export
setClass("GroupedIntegerRanges",
         slot = c(groups = "list",
                  inx = "IntegerList"),
         contains = "DelegatingIntegerRanges",
         validity = validGroupedRanges)

new_grouped_ir <- function(rng, ...) {
    groupings <- make_group_inx(rng, ...)
    # instantiate class
    new("GroupedIntegerRanges",
        delegate = rng, 
        groups = groupings$groups,
        inx = groupings$inx)
}

setMethod("initialize", "GroupedIntegerRanges", 
    function(.Object, delegate, groups, inx, ...) {
        initialize_GroupedRanges(.Object, delegate, groups, inx)
        
})

setMethod("show", "GroupedIntegerRanges", function(object) {
    show_GroupedRanges(object)
})