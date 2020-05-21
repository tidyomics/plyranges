validGroupedRanges <- function(object) {
  group_names <- colnames(object@group_keys)
  check_valid_groups <- !(group_names %in% tbl_vars(object))
  
  if (any(check_valid_groups)) {
    paste("Invalid groups slot:",
          paste(group_names[check_valid_groups], collapse = ","),
          "not found in data.")
  }
}

#' @rdname group_by-ranges
#' @importFrom IRanges extractList 
#' @importFrom S4Vectors splitAsList
#' @importFrom BiocGenerics unlist
#' @export
setClass("GroupedGenomicRanges",
         slot = c(group_keys = "DFrame", 
                  group_indices = "Rle",
                  n = "integer"),
         contains = c("DelegatingGenomicRanges"),
         validity = 
)

initialize_GroupedRanges <- function(.Object, delegate, group_keys, group_indices, n ) {
  .Object@delegate <- delegate
  .Object@group_keys <- group_keys
  .Object@group_indices <- group_indices
  .Object@n <- n
  .Object
}

#' @importFrom methods setMethod initialize
setMethod("initialize", "GroupedGenomicRanges",
          function(.Object, delegate = GRanges(), group_keys = DataFrame(), group_indices = Rle(), n = integer()) {
            initialize_GroupedRanges(.Object, delegate, group_keys, group_indices, n)
          }
)

show_GroupedRanges <- function(object) {
  groups <- colnames(object@group_keys)
  groups <- paste(groups, collapse = ", ")
  output <- c("", utils::capture.output(show(object@delegate)))
  output[1] <- output[2]
  output[2] <- paste("Groups:", groups, paste0("[", object@n, "]"))
  cat(output, sep = "\n")
}
setMethod("show", "GroupedGenomicRanges", function(object) { 
  show_GroupedRanges(object)
})

# --- group-by backend ---
# generates index for grouping variables
new_grouping <- function(rng,  ..., target = "GroupedGenomicRanges") {
  new_groups <- rlang::enquos(...)
  if (length(new_groups) == 0) return(rng)
  
  new_groups <- rlang::quos_auto_name(new_groups)
  # check if we need to mutate, i.e. if quosure is a call
  update_groups <- Filter(rlang::quo_is_call, new_groups)
  
  if (length(update_groups) > 0) {
    rng <- mutate(rng, !!!update_groups)
  }
  
  
  check_names <- !(names(new_groups) %in% tbl_vars(rng))
  
  if (any(check_names)) {
    stop(paste0("Column `", 
               names(new_groups)[check_names],
               "` is unknown"))
  }
  
  group_df <- select(rng, !!!rlang::syms(names(new_groups)), .drop_ranges = TRUE)
  
  unique <- BiocGenerics::unique(group_df)
  inx <- Rle(BiocGenerics::match(group_df, unique))
  n <- nrow(unique)
  new(target, rng, unique, inx, n)
}

# --- GroupedIntegerRanges ---
#' @rdname group_by-ranges
#' @export
setClass("GroupedIntegerRanges",
         slot = c(group_keys = "DFrame", 
                  group_indices = "Rle",
                  n = "integer"),
         contains = "DelegatingIntegerRanges",
         validity = validGroupedRanges)


setMethod("initialize", "GroupedIntegerRanges", 
          function(.Object, delegate = IRanges(), group_keys = DataFrame(), group_indices = Rle(), n = integer()) {
            initialize_GroupedRanges(.Object, delegate, group_keys,group_indices, n)
            
          })

setMethod("show", "GroupedIntegerRanges", function(object) {
  show_GroupedRanges(object)
})



#' @importFrom dplyr tbl_vars
#' @export
tbl_vars.Ranges <- function(x) {
  c("start", "end", "width", names(mcols(x)))
}

#' @export
tbl_vars.GenomicRanges <- function(x) {
  c("start", "end", "width", "strand", "seqnames", names(mcols(x)))
}
