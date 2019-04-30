#' plyranges: a grammar of genomic data manipulation
#'
#' plyranges is a dplyr like API to the Ranges/GenomicRanges infrastructure
#' in Bioconductor.
#'
#' plryanges provides a consistent interface for importing and
#' wrangling genomics data from a variety of sources. The package defines a
#' grammar of genomic data manipulation through a set of verbs. These verbs
#' can be used to construct human readable analysis pipelines based on Ranges
#' objects.
#'
#'   * Modify genomic regions with the `set_width()` and `stretch()` functions.
#'   * Modify genomic regions while fixing the start/end/center coordinates
#'     with the `anchors()` family of functions.
#'   * Sort genomic ranges with `arrange()`.
#'   * Modify, subset, and aggregate genomic data with the `mutate()`,
#'    `filter()`, and `summarise()`functions.
#'   * Any of the above operations can be performed on partitions of the
#'     data with `group_by()`.
#'   * Find nearest neighbour genomic regions with the `join_nearest()` family
#'   of functions.
#'   * Find overlaps between ranges with the `join_overlap_inner()` family of functions.
#'   * Merge all overlapping and adjacent genomic regions with `reduce_ranges()`.
#'   * Merge the end points of all genomic regions with `disjoin_ranges()`.
#'   * Import and write common genomic data formats with the `read_/write_` family
#'    of functions.
#'
#'    For more details on the features of plryanges, read the vignette:
#'    `browseVignettes(package = "plyranges")`
"_PACKAGE"
