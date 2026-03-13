# plyranges: a grammar of genomic data manipulation

plyranges is a dplyr like API to the Ranges/GenomicRanges infrastructure
in Bioconductor.

## Details

plryanges provides a consistent interface for importing and wrangling
genomics data from a variety of sources. The package defines a grammar
of genomic data manipulation through a set of verbs. These verbs can be
used to construct human readable analysis pipelines based on Ranges
objects.

- Modify genomic regions with the
  [`set_width()`](https://tidyomics.github.io/plyranges/reference/ranges-setters.md)
  and
  [`stretch()`](https://tidyomics.github.io/plyranges/reference/stretch.md)
  functions.

- Modify genomic regions while fixing the start/end/center coordinates
  with the `anchors()` family of functions.

- Sort genomic ranges with
  [`arrange()`](https://dplyr.tidyverse.org/reference/arrange.html).

- Modify, subset, and aggregate genomic data with the
  [`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html),
  [`filter()`](https://dplyr.tidyverse.org/reference/filter.html), and
  [`summarise()`](https://dplyr.tidyverse.org/reference/summarise.html)functions.

- Any of the above operations can be performed on partitions of the data
  with
  [`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html).

- Find nearest neighbour genomic regions with the
  [`join_nearest()`](https://tidyomics.github.io/plyranges/reference/ranges-nearest.md)
  family of functions.

- Find overlaps between ranges with the
  [`join_overlap_inner()`](https://tidyomics.github.io/plyranges/reference/overlap-joins.md)
  family of functions.

- Merge all overlapping and adjacent genomic regions with
  [`reduce_ranges()`](https://tidyomics.github.io/plyranges/reference/ranges-reduce.md).

- Merge the end points of all genomic regions with
  [`disjoin_ranges()`](https://tidyomics.github.io/plyranges/reference/ranges-disjoin.md).

- Import and write common genomic data formats with the `read_/write_`
  family of functions.

For more details on the features of plryanges, read the vignette:
`browseVignettes(package = "plyranges")`

## See also

Useful links:

- Report bugs at <https://github.com/tidyomics/plyranges>

## Author

**Maintainer**: Michael Love <michaelisaiahlove@gmail.com>
\[contributor\]

Authors:

- Stuart Lee ([ORCID](https://orcid.org/0000-0003-1179-8436))

- Michael Lawrence \[contributor\]

- Dianne Cook \[contributor\]

Other contributors:

- Spencer Nystrom ([ORCID](https://orcid.org/0000-0003-1000-1579))
  \[contributor\]

- Pierre-Paul Axisa \[contributor\]
