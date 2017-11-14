# Introduction to the plyranges package
{:.no_toc}

<p class="author-name">Stuart Lee</p>

## Contents
{:.no_toc}

-   [Introduction](#introduction)
-   [Ranges revisted](#ranges-revisted)
-   [Constructing Ranges](#constructing-ranges)
-   [Arithmetic on Ranges](#arithmetic-on-ranges)
-   [Grouping Ranges](#grouping-ranges)
-   [Restricting Ranges](#restricting-ranges)
-   [Joins, or another way at looking at overlaps between Ranges](#joins-or-another-way-at-looking-at-overlaps-between-ranges)
-   [Summarising Ranges](#summarising-ranges)
-   [Data Import/Output](#data-importoutput)
-   [Appendix](#appendix)
{:toc}

## Introduction

The `plyranges` package provides a [fluent](https://en.wikipedia.org/wiki/Fluent_interface) interface for performing common genomic data analysis tasks. One goal of this package is to decrease the learning curve for Bioconductor and S4 classes (especially for new users) by providing a consistent interface to common classes, `IRanges` and `GRanges` classes.

This document provides an overview of some of the features of the `plyranges` grammar for analysing genomic data. I give a brief introduction to the `Ranges` and `GenomicRanges` classes and then detail each feature of the `plyranges` library: construction, arithmetic, restriction, summarisation, joins, and data import/output.

## Ranges revisted

`Ranges` objects can either represent intervals of integers as `IRanges` (which have start, end and width attributes) or represent genomic intervals (which have additional attributes, sequence name, and strand) as `GRanges`. In addition, both types of `Ranges` can store information about their intervals as metadata columns (for example GC content over a genomic interval).

`Ranges` objects follow the tidy data principle: each row of a `Ranges` object corresponds to an interval, while each column will represent a variable about that interval, and generally each object will represent a single unit of observation (like gene annotations).

Consequently, `Ranges` objects provide a powerful representation for reasoning about genomic data. In this vignette, you will learn more about `Ranges` objects and how via grouping, restriction and summarisation you can perform common data tasks.

## Constructing Ranges

To construct an `IRanges` we require there are at least two columns that represent at either a starting coordinate, finishing coordinate or the width of the interval. To construct a `GRanges` we require a column that represents that sequence name ( contig or chromosome id), and an optional column to represent the strandedness of an interval.

``` r
suppressPackageStartupMessages(library(plyranges))
set.seed(100)
df <- data.frame(start=c(2:-1, 13:15), 
                 width=c(0:3, 2:0))

## produces IRanges
rng <- df %>% Ranges()
rng
```

    ## IRanges object with 7 ranges and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]         2         1         0
    ##   [2]         1         1         1
    ##   [3]         0         1         2
    ##   [4]        -1         1         3
    ##   [5]        13        14         2
    ##   [6]        14        14         1
    ##   [7]        15        14         0

``` r

## seqname is required for GRanges, metadata is automatically kept
grng <- df %>% 
  transform(seqnames = sample(c("chr1", "chr2"), 7, replace = TRUE),
         strand = sample(c("+", "-"), 7, replace = TRUE),
         gc = runif(7)) %>% 
  Ranges()

grng
```

    ## GRanges object with 7 ranges and 1 metadata column:
    ##       seqnames    ranges strand |                gc
    ##          <Rle> <IRanges>  <Rle> |         <numeric>
    ##   [1]     chr1  [ 2,  1]      + |  0.76255108229816
    ##   [2]     chr1  [ 1,  1]      - | 0.669021712383255
    ##   [3]     chr2  [ 0,  1]      + | 0.204612161964178
    ##   [4]     chr1  [-1,  1]      - | 0.357524853432551
    ##   [5]     chr1  [13, 14]      - | 0.359475114848465
    ##   [6]     chr1  [14, 14]      + | 0.690290528349578
    ##   [7]     chr2  [15, 14]      + | 0.535811153938994
    ##   -------
    ##   seqinfo: 2 sequences from an unspecified genome; no seqlengths

## Arithmetic on Ranges

Sometimes you want to modify a genomic interval by altering the width of the interval while leaving the start, end or midpoint of the coordinates unaltered. This is achieved with the `set_width` verb along with `anchor_*` adverbs.

``` r
rng <- Ranges(data.frame(start=c(1, 2, 3), end=c(5, 2, 8)))
grng <- Ranges(data.frame(start=c(1, 2, 3), end=c(5, 2, 8), 
                          seqnames = "seq1",
                          strand = c("+", "*", "-")))
set_width(rng, 10)
```

    ## IRanges object with 3 ranges and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]         1        10        10
    ##   [2]         2        11        10
    ##   [3]         3        12        10

``` r
set_width(anchor_start(rng), 10)
```

    ## IRanges object with 3 ranges and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]         1        10        10
    ##   [2]         2        11        10
    ##   [3]         3        12        10

``` r
set_width(anchor_end(rng), 10)
```

    ## IRanges object with 3 ranges and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]        -4         5        10
    ##   [2]        -7         2        10
    ##   [3]        -1         8        10

``` r
set_width(anchor_center(rng), 10)
```

    ## IRanges object with 3 ranges and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]        -2         7        10
    ##   [2]        -3         6        10
    ##   [3]         1        10        10

``` r
set_width(anchor_3p(grng), 10) # leave negative strand fixed
```

    ## GRanges object with 3 ranges and 0 metadata columns:
    ##       seqnames    ranges strand
    ##          <Rle> <IRanges>  <Rle>
    ##   [1]     seq1   [1, 10]      +
    ##   [2]     seq1   [2,  2]      *
    ##   [3]     seq1   [3,  8]      -
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome; no seqlengths

``` r
set_width(anchor_5p(grng), 10) # leave positve strand fixed
```

    ## GRanges object with 3 ranges and 0 metadata columns:
    ##       seqnames    ranges strand
    ##          <Rle> <IRanges>  <Rle>
    ##   [1]     seq1   [ 1, 5]      +
    ##   [2]     seq1   [ 2, 2]      *
    ##   [3]     seq1   [-1, 8]      -
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome; no seqlengths

Similarly, you can modify the width of an interval using the `stretch` verb. Without anchoring, this function will extend the interval in either direction by an integer amount. With anchoring, either the start, end or midpoint are preserved.

``` r
rng2 <- stretch(anchor_center(rng), 10)
rng2
```

    ## IRanges object with 3 ranges and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]        -7        13        21
    ##   [2]        -8        12        21
    ##   [3]        -5        16        22

``` r
stretch(anchor_end(rng2), 10)
```

    ## IRanges object with 3 ranges and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]         3        13        11
    ##   [2]         2        12        11
    ##   [3]         5        16        12

``` r
stretch(anchor_start(rng2), 10)
```

    ## IRanges object with 3 ranges and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]        -7        23        31
    ##   [2]        -8        22        31
    ##   [3]        -5        26        32

``` r
stretch(anchor_3p(grng), 10)
```

    ## GRanges object with 3 ranges and 0 metadata columns:
    ##       seqnames    ranges strand
    ##          <Rle> <IRanges>  <Rle>
    ##   [1]     seq1  [-9, 15]      +
    ##   [2]     seq1  [ 2,  2]      *
    ##   [3]     seq1  [ 3,  8]      -
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome; no seqlengths

``` r
stretch(anchor_5p(grng), 10)
```

    ## GRanges object with 3 ranges and 0 metadata columns:
    ##       seqnames    ranges strand
    ##          <Rle> <IRanges>  <Rle>
    ##   [1]     seq1  [ 1,  5]      +
    ##   [2]     seq1  [ 2,  2]      *
    ##   [3]     seq1  [-7, 18]      -
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome; no seqlengths

`*Ranges` can be shifted left or right. If strand information is available we can also shift upstream or downstream.

``` r
shift_left(rng, 10)
```

    ## IRanges object with 3 ranges and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]        -9        -5         5
    ##   [2]        -8        -8         1
    ##   [3]        -7        -2         6

``` r
shift_right(rng, 10)
```

    ## IRanges object with 3 ranges and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]        11        15         5
    ##   [2]        12        12         1
    ##   [3]        13        18         6

``` r
shift_upstream(grng, 10)
```

    ## GRanges object with 3 ranges and 0 metadata columns:
    ##       seqnames    ranges strand
    ##          <Rle> <IRanges>  <Rle>
    ##   [1]     seq1  [-9, -5]      +
    ##   [2]     seq1  [ 2,  2]      *
    ##   [3]     seq1  [13, 18]      -
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome; no seqlengths

``` r
shift_downstream(grng, 10)
```

    ## GRanges object with 3 ranges and 0 metadata columns:
    ##       seqnames    ranges strand
    ##          <Rle> <IRanges>  <Rle>
    ##   [1]     seq1  [11, 15]      +
    ##   [2]     seq1  [ 2,  2]      *
    ##   [3]     seq1  [-7, -2]      -
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome; no seqlengths

## Grouping Ranges

`plyranges` introduces a new class of `Ranges` called `RangesGrouped`, this is a similiar idea to the grouped `data.frame\tibble` in `dplyr`. Note that grouping does not change the structure of a `Ranges` object only how certain functions act on it.

Grouping can act on either the core components or the metadata columns.

``` r
grng <- data.frame(seqnames = sample(c("chr1", "chr2"), 7, replace = TRUE),
         strand = sample(c("+", "-"), 7, replace = TRUE),
         gc = runif(7),
         start = 1:7,
         width = 10) %>%
  Ranges()

grng_by_strand <- grng %>%
  group_by(strand)

grng_by_strand
```

    ## GRangesGrouped object with 7 ranges and 1 metadata column:
    ## Groups: strand
    ##       seqnames    ranges strand |                gc
    ##          <Rle> <IRanges>  <Rle> |         <numeric>
    ##   [1]     chr2   [1, 10]      - | 0.889453538926318
    ##   [2]     chr2   [2, 11]      + | 0.180407245177776
    ##   [3]     chr2   [3, 12]      + | 0.629390850430354
    ##   [4]     chr1   [4, 13]      - | 0.989564136601985
    ##   [5]     chr1   [5, 14]      + | 0.130288870073855
    ##   [6]     chr2   [6, 15]      - | 0.330660525709391
    ##   [7]     chr2   [7, 16]      - | 0.865120546659455
    ##   -------
    ##   seqinfo: 2 sequences from an unspecified genome; no seqlengths

## Restricting Ranges

The verb `filter` can be used to restrict rows in the Ranges. Note that grouping will cause the `filter` to act within each group of the data.

``` r
grng %>% filter(gc < 0.3)
```

    ## GRanges object with 2 ranges and 1 metadata column:
    ##       seqnames    ranges strand |                gc
    ##          <Rle> <IRanges>  <Rle> |         <numeric>
    ##   [1]     chr2   [2, 11]      + | 0.180407245177776
    ##   [2]     chr1   [5, 14]      + | 0.130288870073855
    ##   -------
    ##   seqinfo: 2 sequences from an unspecified genome; no seqlengths

``` r

## filtering by group 
grng_by_strand %>% filter(gc == max(gc)) 
```

    ## GRangesGrouped object with 2 ranges and 1 metadata column:
    ## Groups: 
    ##       seqnames    ranges strand |                gc
    ##          <Rle> <IRanges>  <Rle> |         <numeric>
    ##   [1]     chr2   [3, 12]      + | 0.629390850430354
    ##   [2]     chr1   [4, 13]      - | 0.989564136601985
    ##   -------
    ##   seqinfo: 2 sequences from an unspecified genome; no seqlengths

We also provide the convience methods `filter_by_overlaps` and `filter_by_non_overlaps` for restricting by overlapping ranges.

``` r
ir0 <- data.frame(start = c(5,10, 15,20), width = 5) %>%
  Ranges()
ir1 <- data.frame(start = 2:6, width = 3:7) %>%
  Ranges()
ir0
```

    ## IRanges object with 4 ranges and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]         5         9         5
    ##   [2]        10        14         5
    ##   [3]        15        19         5
    ##   [4]        20        24         5

``` r
ir1
```

    ## IRanges object with 5 ranges and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]         2         4         3
    ##   [2]         3         6         4
    ##   [3]         4         8         5
    ##   [4]         5        10         6
    ##   [5]         6        12         7

``` r
ir0 %>% filter_by_overlaps(ir1)
```

    ## IRanges object with 2 ranges and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]         5         9         5
    ##   [2]        10        14         5

``` r
ir0 %>% filter_by_non_overlaps(ir1) 
```

    ## IRanges object with 2 ranges and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]        15        19         5
    ##   [2]        20        24         5

## Joins, or another way at looking at overlaps between Ranges

We can think about finding overlaps as a type of join between two `*Ranges` where we are using the intervals as keys. Of most use are the `join_overlap_inner` and `join_overlap_left` functions. The former will return a `Ranges` object where the start, end, and width coordinates correspond to the amount of any overlap. By default, the width of the input ranges is added as metadata. This allows downstreaming filtering by the amount of overlap.

``` r
join_overlap_inner(ir0, ir1)
```

    ## IRanges object with 6 ranges and 2 metadata columns:
    ##           start       end     width |   width.x   width.y
    ##       <integer> <integer> <integer> | <integer> <integer>
    ##   [1]         5         6         2 |         5         4
    ##   [2]         5         8         4 |         5         5
    ##   [3]         5         9         5 |         5         6
    ##   [4]         6         9         4 |         5         7
    ##   [5]        10        10         1 |         5         6
    ##   [6]        10        12         3 |         5         7

The `join_overlap_left` method will add the coordinates of the right hand ranges that overlap those ranges on the left. Compared with `filter_by_overlaps` above, the overlap join expands the ranges to give information about each range on the right hand side that overlaps those on the left. We also provide a convienence method called `find_overlaps` that computes the same result.

``` r
join_overlap_left(ir0, ir1)
```

    ## IRanges object with 6 ranges and 3 metadata columns:
    ##           start       end     width |   start.y     end.y   width.y
    ##       <integer> <integer> <integer> | <integer> <integer> <integer>
    ##   [1]         5         9         5 |         3         6         4
    ##   [2]         5         9         5 |         4         8         5
    ##   [3]         5         9         5 |         5        10         6
    ##   [4]         5         9         5 |         6        12         7
    ##   [5]        10        14         5 |         5        10         6
    ##   [6]        10        14         5 |         6        12         7

We also provide methods for finding, nearest, preceding or following ranges. These methods nest a ranges column upon returning (maybe better to return just coordinates?).

``` r
join_nearest(ir0, ir1)
```

    ## IRanges object with 4 ranges and 1 metadata column:
    ##           start       end     width |   nearest
    ##       <integer> <integer> <integer> | <IRanges>
    ##   [1]         5         9         5 |   [6, 12]
    ##   [2]        10        14         5 |   [6, 12]
    ##   [3]        15        19         5 |   [6, 12]
    ##   [4]        20        24         5 |   [6, 12]

``` r
join_follow(ir0, ir1)
```

    ## IRanges object with 4 ranges and 1 metadata column:
    ##           start       end     width |   follows
    ##       <integer> <integer> <integer> | <IRanges>
    ##   [1]         5         9         5 |   [2,  4]
    ##   [2]        10        14         5 |   [4,  8]
    ##   [3]        15        19         5 |   [6, 12]
    ##   [4]        20        24         5 |   [6, 12]

``` r
join_precede(ir0, ir1) # nothing precedes returns empty ranges
```

    ## IRanges object with 0 ranges and 1 metadata column:
    ##        start       end     width |  precedes
    ##    <integer> <integer> <integer> | <IRanges>

``` r
join_precede(ir1, ir0)
```

    ## IRanges object with 5 ranges and 1 metadata column:
    ##           start       end     width |  precedes
    ##       <integer> <integer> <integer> | <IRanges>
    ##   [1]         2         4         3 |  [ 5,  9]
    ##   [2]         3         6         4 |  [10, 14]
    ##   [3]         4         8         5 |  [10, 14]
    ##   [4]         5        10         6 |  [15, 19]
    ##   [5]         6        12         7 |  [15, 19]

For `GRanges` objects by default strand is not considered when performing overlap joins. To include strand when finding overlaps use append the prefix `directed`. To restrict overlapping ranges to those within the query range use the prefix `within`.

It's also possible to group by overlaps. Using this approach we can count the number of overlaps that are greater than 0.

``` r
grp_by_olap <- ir0 %>% 
  group_by_overlaps(ir1)
grp_by_olap
```

    ## IRangesGrouped object with 6 ranges and 4 metadata columns:
    ## Groups: query
    ##           start       end     width | start.subject end.subject width.subject
    ##       <integer> <integer> <integer> |     <integer>   <integer>     <integer>
    ##   [1]         5         9         5 |             3           6             4
    ##   [2]         5         9         5 |             4           8             5
    ##   [3]         5         9         5 |             5          10             6
    ##   [4]         5         9         5 |             6          12             7
    ##   [5]        10        14         5 |             5          10             6
    ##   [6]        10        14         5 |             6          12             7
    ##           query
    ##       <integer>
    ##   [1]         1
    ##   [2]         1
    ##   [3]         1
    ##   [4]         1
    ##   [5]         2
    ##   [6]         2

``` r
grp_by_olap %>%
  mutate(n_overlaps = length(query))
```

    ## IRangesGrouped object with 6 ranges and 5 metadata columns:
    ## Groups: query
    ##           start       end     width | start.subject end.subject width.subject
    ##       <integer> <integer> <integer> |     <integer>   <integer>     <integer>
    ##   [1]         5         9         5 |             3           6             4
    ##   [2]         5         9         5 |             4           8             5
    ##   [3]         5         9         5 |             5          10             6
    ##   [4]         5         9         5 |             6          12             7
    ##   [5]        10        14         5 |             5          10             6
    ##   [6]        10        14         5 |             6          12             7
    ##           query n_overlaps
    ##       <integer>  <integer>
    ##   [1]         1          4
    ##   [2]         1          4
    ##   [3]         1          4
    ##   [4]         1          4
    ##   [5]         2          2
    ##   [6]         2          2

Of course we can also add overlap counts via the `count_overlaps` function.

``` r
ir0 %>%
  mutate(n_overlaps = count_overlaps(., ir1))
```

    ## IRanges object with 4 ranges and 1 metadata column:
    ##           start       end     width | n_overlaps
    ##       <integer> <integer> <integer> |  <integer>
    ##   [1]         5         9         5 |          4
    ##   [2]        10        14         5 |          2
    ##   [3]        15        19         5 |          0
    ##   [4]        20        24         5 |          0

## Summarising Ranges

The `summarise` function will return a tibble because the information required to return a \_\*Ranges\_ object is lost. It is often most useful to use `summarise` in combination with the `group_by` family of functions.

``` r
ir1 <- ir1 %>%
  mutate(gc = runif(length(.)))

ir0 %>% 
  group_by_overlaps(ir1) %>%
  summarise(gc = mean(gc))
```

    ## # A tibble: 2 x 2
    ##   query        gc
    ##   <int>     <dbl>
    ## 1     1 0.6755545
    ## 2     2 0.6357952

## Data Import/Output

We provide convienence functions via `rtracklayer` for reading/writing the following data formats from/to \_\*Ranges\_ objects.

-   BED: `read/write_bed`
-   BEDGraph: `read/write_bedgraph`
-   GFF(1-3): `read/write_gff(1-3)`
-   BigWig: `read/write_bw`
-   Wig: `read/write_wig`
-   narrowPeaks: `read/write_narrowpeaks`

## Appendix

``` r
sessionInfo()
```

    ## R version 3.4.2 (2017-09-28)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Sierra 10.12.6
    ## 
    ## Matrix products: default
    ## BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ## [1] plyranges_0.1.0      BiocStyle_2.6.0      GenomicRanges_1.30.0
    ## [4] GenomeInfoDb_1.14.0  IRanges_2.12.0       S4Vectors_0.16.0    
    ## [7] BiocGenerics_0.24.0 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] SummarizedExperiment_1.8.0 tidyselect_0.2.3          
    ##  [3] purrr_0.2.4                lattice_0.20-35           
    ##  [5] testthat_1.0.2             htmltools_0.3.6           
    ##  [7] rtracklayer_1.38.0         yaml_2.1.14               
    ##  [9] XML_3.98-1.9               rlang_0.1.4               
    ## [11] glue_1.2.0                 withr_2.1.0               
    ## [13] BiocParallel_1.12.0        bindrcpp_0.2              
    ## [15] matrixStats_0.52.2         GenomeInfoDbData_0.99.1   
    ## [17] bindr_0.1                  stringr_1.2.0             
    ## [19] zlibbioc_1.24.0            Biostrings_2.46.0         
    ## [21] commonmark_1.4             devtools_1.13.4           
    ## [23] memoise_1.1.0              evaluate_0.10.1           
    ## [25] Biobase_2.38.0             knitr_1.17                
    ## [27] BiocInstaller_1.28.0       Rcpp_0.12.13              
    ## [29] backports_1.1.1            DelayedArray_0.4.1        
    ## [31] desc_1.1.1                 XVector_0.18.0            
    ## [33] Rsamtools_1.30.0           digest_0.6.12             
    ## [35] stringi_1.1.5              bookdown_0.5              
    ## [37] dplyr_0.7.4                grid_3.4.2                
    ## [39] rprojroot_1.2              tools_3.4.2               
    ## [41] bitops_1.0-6               magrittr_1.5              
    ## [43] RCurl_1.95-4.8             tibble_1.3.4              
    ## [45] crayon_1.3.4               tidyr_0.7.2               
    ## [47] pkgconfig_2.0.1            Matrix_1.2-11             
    ## [49] xml2_1.1.1                 assertthat_0.2.0          
    ## [51] rmarkdown_1.7              roxygen2_6.0.1            
    ## [53] rstudioapi_0.7             R6_2.2.2                  
    ## [55] GenomicAlignments_1.14.0   compiler_3.4.2
