<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build Status](https://travis-ci.org/sa-lee/plyranges.svg?branch=master)](https://travis-ci.org/sa-lee/plyranges)

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

### Modifying Genomic Intervals

Sometimes you want to modify a genomic interval by altering it's coordinateswidth and leave the start, end or midpoint of the coordinates unaltered. This is achieved with the `set_width` verb along with `anchor_*` adverbs.

``` r
rng <- Ranges(df)
set_width(rng, 100L)
#> IRanges object with 7 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]         2       101       100
#>   [2]         1       100       100
#>   [3]         0        99       100
#>   [4]        -1        98       100
#>   [5]        13       112       100
#>   [6]        14       113       100
#>   [7]        15       114       100
set_width(anchor_start(rng), 100L)
#> IRanges object with 7 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]         2       101       100
#>   [2]         1       100       100
#>   [3]         0        99       100
#>   [4]        -1        98       100
#>   [5]        13       112       100
#>   [6]        14       113       100
#>   [7]        15       114       100
set_width(anchor_end(rng), 100L)
#> IRanges object with 7 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]       -98         1       100
#>   [2]       -98         1       100
#>   [3]       -98         1       100
#>   [4]       -98         1       100
#>   [5]       -85        14       100
#>   [6]       -85        14       100
#>   [7]       -85        14       100
set_width(anchor_center(rng), 100L)
#> IRanges object with 7 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]       -48        51       100
#>   [2]       -49        50       100
#>   [3]       -49        50       100
#>   [4]       -50        49       100
#>   [5]       -36        63       100
#>   [6]       -36        63       100
#>   [7]       -35        64       100
```

Similarly, you can modify the width of an interval using the `stretch` verb. Without anchoring, this function will extend the interval in either direction by an integer amount. With anchoring, either the start, end or midpoint are preserved.

``` r
rng2 <- stretch(anchor_center(rng), 10)
rng2
#> IRanges object with 7 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]        -9        12        22
#>   [2]        -9        11        21
#>   [3]       -10        11        22
#>   [4]       -10        10        21
#>   [5]         3        24        22
#>   [6]         4        24        21
#>   [7]         4        25        22
stretch(anchor_end(rng2), 1L)
#> IRanges object with 7 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]        -8        12        21
#>   [2]        -8        11        20
#>   [3]        -9        11        21
#>   [4]        -9        10        20
#>   [5]         4        24        21
#>   [6]         5        24        20
#>   [7]         5        25        21
stretch(anchor_start(rng2), 10L)
#> IRanges object with 7 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]        -9        22        32
#>   [2]        -9        21        31
#>   [3]       -10        21        32
#>   [4]       -10        20        31
#>   [5]         3        34        32
#>   [6]         4        34        31
#>   [7]         4        35        32
```

*Ranges* can be shifted left or right.

``` r
shift_left(rng, 10L)
#> IRanges object with 7 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]        -8        -9         0
#>   [2]        -9        -9         1
#>   [3]       -10        -9         2
#>   [4]       -11        -9         3
#>   [5]         3         4         2
#>   [6]         4         4         1
#>   [7]         5         4         0
shift_right(rng, 10L)
#> IRanges object with 7 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]        12        11         0
#>   [2]        11        11         1
#>   [3]        10        11         2
#>   [4]         9        11         3
#>   [5]        23        24         2
#>   [6]        24        24         1
#>   [7]        25        24         0
```

[1] tentative title
