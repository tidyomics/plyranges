# plyranges: fluent genomic data analysis

*plyranges* provides a consistent interface for importing and wrangling
genomics data from a variety of sources. The package defines a grammar
of genomic data transformation based on *dplyr* and the Bioconductor
packages *IRanges*, *GenomicRanges*, and *rtracklayer*. It does this by
providing a set of verbs for developing analysis pipelines based on
*GRanges* objects that represent genomic regions:

- Modify genomic regions with the `mutate()` and
  [`stretch()`](https://tidyomics.github.io/plyranges/reference/stretch.md)
  functions.
- Modify genomic regions while fixing the start/end/center coordinates
  with the `anchor_` family of functions.
- Sort genomic ranges with `arrange()`.
- Modify, subset, and aggregate genomic data with the `mutate()`,
  [`filter()`](https://rdrr.io/r/stats/filter.html), and
  `summarise()`functions.
- Any of the above operations can be performed on partitions of the data
  with `group_by()`.
- Find nearest neighbour genomic regions with the `join_nearest_` family
  of functions.
- Find overlaps between ranges with the `join_overlaps_` family of
  functions.
- Add additional metadata between ranges and a table with the
  `join_mcols_` family of functions.
- Merge all overlapping and adjacent genomic regions with
  [`reduce_ranges()`](https://tidyomics.github.io/plyranges/reference/ranges-reduce.md).
- Merge the end points of all genomic regions with
  [`disjoin_ranges()`](https://tidyomics.github.io/plyranges/reference/ranges-disjoin.md).
- Import and write common genomic data formats with the `read_/write_`
  family of functions.

# Documentation

For more details on the features of *plyranges*, read the [introductory
vignette](https://tidyomics.github.io/plyranges/articles/an-introduction.html)
and the [examples
vignette](https://tidyomics.github.io/plyranges/articles/more-examples.html)

For a complete case-study on using *plyranges* to combine ATAC-seq and
RNA-seq results read the [*fluentGenomics*
workflow](https://tidyomics.github.io/fluentGenomics).

*plyranges* is part of the [tidyomics](https://github.com/tidyomics)
project, providing a `dplyr`-based interface for many types of genomics
datasets represented in Bioconductor.

# Installation

*plyranges* can be installed from the latest Bioconductor release:

``` r

# install.packages("BiocManager")
BiocManager::install("plyranges")
```

To install the development version from GitHub:

``` r

BiocManager::install("tidyomics/plyranges")
```

# Learning more

In addition to the two package vignettes, see the following for more
informtion:

- The [fluentGenomics workflow](https://sa-lee.github.io/fluentGenomics)
  package shows how to combine differential gene expression and
  differential chromatin accessibility using *plyranges*.

- The [extended vignette in the plyrangesWorkshops
  package](https://github.com/sa-lee/plyrangesWorkshops) has a detailed
  walk through of using *plyranges* for coverage analysis.

- The collection of genomic range applications including *plyranges*:
  [tidy ranges
  tutorial](https://tidyomics.github.io/tidy-ranges-tutorial).

# Citation

If you found *plyranges* useful for your work please cite our
[paper](http://doi.org/10.1186/s13059-018-1597-8):

``` R
@ARTICLE{Lee2019,
  title    = "plyranges: a grammar of genomic data transformation",
  author   = "Lee, Stuart and Cook, Dianne and Lawrence, Michael",
  journal  = "Genome Biol.",
  volume   =  20,
  number   =  1,
  pages    = "4",
  month    =  jan,
  year     =  2019,
  url      = "http://dx.doi.org/10.1186/s13059-018-1597-8",
  doi      = "10.1186/s13059-018-1597-8",
  pmc      = "PMC6320618"
}
```

# Contributing

We welcome contributions from the R/Bioconductor community. We ask that
contributors follow the [code of
conduct](https://tidyomics.github.io/plyranges/CODE_OF_CONDUCT.md) and
the guide outlined
[here](https://tidyomics.github.io/plyranges/CONTRIBUTING.md).
