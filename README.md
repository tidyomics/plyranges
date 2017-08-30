<!-- README.md is generated from README.Rmd. Please edit that file -->
plyranges
=========

The `plyranges` [1] package provides a [fluent](https://en.wikipedia.org/wiki/Fluent_interface) interface for performing common genomic data analysis tasks. One goal of this package is to decrease the learning curve for Bioconductor and S4 classes (especially for new users) by providing a consistent interface to common classes, *IRanges*, *GRanges* and *SummarizedExperiment*.

To achieve this we have constructed a grammar of genomic data manipulation based on `dplyr` and the Bioconductor packages that define the aformentioned classes. All the methods defined in `plyranges` are type-consistent and chainable. As a consequence, users familiar with `dplyr` should be able to read and understand `plyranges` code without difficulty.

*Note this package is under heavy development*

Installation
============

Currently only the development version is available

``` r
source("https://bioconductor.org/biocLite.R")
biocLite("sa-lee/plyranges")
```

Constructing Ranges
===================

Ranges can be easily constructed from a data.frame or tibble. To construct an `IRanges` we require there are at least two columns that represent at either a starting coordinate, finishing coordinate or the width of the interval. To construct a `GRanges` we require a column that represents that sequence name ( contig or chromosome id), and an optional column to represent the strandedness of an interval.

``` r
set.seed(100)
df <- data.frame(start=c(2:-1, 13:15), 
                 width=c(0:3, 2:0))

# produces IRanges
df %>% Ranges()
#> IRanges object with 7 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]         2         1         0
#>   [2]         1         1         1
#>   [3]         0         1         2
#>   [4]        -1         1         3
#>   [5]        13        14         2
#>   [6]        14        14         1
#>   [7]        15        14         0

# produces GRanges
df %>% 
  mutate(seqnames = sample(c("chr1", "chr2"), 7, replace = TRUE)) %>% 
  Ranges()
#> GRanges object with 7 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1  [ 2,  1]      *
#>   [2]     chr1  [ 1,  1]      *
#>   [3]     chr2  [ 0,  1]      *
#>   [4]     chr1  [-1,  1]      *
#>   [5]     chr1  [13, 14]      *
#>   [6]     chr1  [14, 14]      *
#>   [7]     chr2  [15, 14]      *
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

df  %>% 
  mutate(strand = sample(c("+", "-"), 7, replace = TRUE)) %>% 
  Ranges()
#> IRanges object with 7 ranges and 1 metadata column:
#>           start       end     width |      strand
#>       <integer> <integer> <integer> | <character>
#>   [1]         2         1         0 |           +
#>   [2]         1         1         1 |           -
#>   [3]         0         1         2 |           +
#>   [4]        -1         1         3 |           -
#>   [5]        13        14         2 |           -
#>   [6]        14        14         1 |           +
#>   [7]        15        14         0 |           +

# seqname is required for GRanges, metadata is automatically kept
df %>% 
  mutate(seqnames = sample(c("chr1", "chr2"), 7, replace = TRUE),
         strand = sample(c("+", "-"), 7, replace = TRUE),
         gc = runif(7)) %>% 
  Ranges()
#> GRanges object with 7 ranges and 1 metadata column:
#>       seqnames    ranges strand |                gc
#>          <Rle> <IRanges>  <Rle> |         <numeric>
#>   [1]     chr2  [ 2,  1]      - | 0.549096710281447
#>   [2]     chr2  [ 1,  1]      - | 0.277723756618798
#>   [3]     chr1  [ 0,  1]      - | 0.488305994076654
#>   [4]     chr1  [-1,  1]      + | 0.928505074931309
#>   [5]     chr1  [13, 14]      + | 0.348691981751472
#>   [6]     chr2  [14, 14]      - | 0.954157707514241
#>   [7]     chr2  [15, 14]      - | 0.695274139055982
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```

[1] tentative title
