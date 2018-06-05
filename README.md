
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Build
Status](https://travis-ci.org/sa-lee/plyranges.svg?branch=master)](https://travis-ci.org/sa-lee/plyranges)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/sa-lee/plyranges?branch=master&svg=true)](https://ci.appveyor.com/project/sa-lee/plyranges)
[![Coverage
Status](https://img.shields.io/codecov/c/github/sa-lee/plyranges/master.svg)](https://codecov.io/github/sa-lee/plyranges?branch=master)

# plyranges

`plryanges` provides a consistent interface for importing and wrangling
genomics data from a variety of sources. The package defines a grammar
of genomic data manipulation based on `dplyr` and the Bioconductor
packages `IRanges`, `GenomicRanges`, and `rtracklayer`. It does this by
providing a set of verbs for developing analysis pipelines based on
*Ranges* objects that represent genomic regions:

  - Modify genomic regions with the `mutate()` and `stretch()`
    functions.
  - Modify genomic regions while fixing the start/end/center coordinates
    with the `anchor_` family of functions.
  - Sort genomic ranges with `arrange()`.
  - Modify, subset, and aggregate genomic data with the `mutate()`,
    `filter()`, and `summarise()`functions.
  - Any of the above operations can be performed on partitions of the
    data with `group_by()`.
  - Find nearest neighbour genomic regions with the `join_nearest_`
    family of functions.
  - Find overlaps between ranges with the `join_overlaps_` family of
    functions.
  - Merge all overlapping and adjacent genomic regions with
    `reduce_ranges()`.
  - Merge the end points of all genomic regions with `disjoin_ranges()`.
  - Import and write common genomic data formats with the `read_/write_`
    family of functions.

For more details on the features of plryanges, read the [introduction
vignette](https://sa-lee.github.io/plyranges/articles/an-introduction.html).

# Installation

The package is currently available from Bioconductor.

``` r
source("https://bioconductor.org/biocLite.R")
biocLite("plyranges")
```

To install the development version from GitHub:

``` r
source("https://bioconductor.org/biocLite.R")
biocLite("sa-lee/plyranges")
```

# Quick overview

## About `Ranges`

`Ranges` objects can either represent sets of integers as `IRanges`
(which have start, end and width attributes) or represent genomic
intervals (which have additional attributes, sequence name, and strand)
as `GRanges`. In addition, both types of `Ranges` can store information
about their intervals as metadata columns (for example GC content over a
genomic interval).

`Ranges` objects follow the tidy data principle: each row of a `Ranges`
object corresponds to an interval, while each column will represent a
variable about that interval, and generally each object will represent a
single unit of observation (like gene annotations).

We can construct a `IRanges` object from a `data.frame` with a `start`
or `width` using the `as_iranges()` method.

``` r
library(plyranges)
df <- data.frame(start = 1:5, width = 5)
as_iranges(df)
#> IRanges object with 5 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]         1         5         5
#>   [2]         2         6         5
#>   [3]         3         7         5
#>   [4]         4         8         5
#>   [5]         5         9         5
# alternatively with end
df <- data.frame(start = 1:5, end = 5:9)
as_iranges(df)
#> IRanges object with 5 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]         1         5         5
#>   [2]         2         6         5
#>   [3]         3         7         5
#>   [4]         4         8         5
#>   [5]         5         9         5
```

We can also construct a `GRanges` object in a similar manner. Note that
a `GRanges` object requires at least a seqnames column to be present in
the data.frame (but not necessarily a strand column).

``` r
df <- data.frame(seqnames = c("chr1", "chr2", "chr2", "chr1", "chr2"),
                 start = 1:5,
                 width = 5)
as_granges(df)
#> GRanges object with 5 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1       1-5      *
#>   [2]     chr2       2-6      *
#>   [3]     chr2       3-7      *
#>   [4]     chr1       4-8      *
#>   [5]     chr2       5-9      *
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
# strand can be specified with `+`, `*` (mising) and `-`
df$strand <- c("+", "+", "-", "-", "*")
as_granges(df)
#> GRanges object with 5 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1       1-5      +
#>   [2]     chr2       2-6      +
#>   [3]     chr2       3-7      -
#>   [4]     chr1       4-8      -
#>   [5]     chr2       5-9      *
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```

# Example: finding GWAS hits that overlap known exons

Let’s look at a more a realistic example (taken from HelloRanges
vignette).

Suppose we have two *GRanges* objects: one containing coordinates of
known exons and another containing SNPs from a GWAS.

The first and last 5 exons are printed below, there are two additional
columns corresponding to the exon name, and a score.

We could check the number of exons per chromosome using `group_by` and
`summarise`.

``` r
exons 
#> GRanges object with 459752 ranges and 2 metadata columns:
#>            seqnames            ranges strand |
#>               <Rle>         <IRanges>  <Rle> |
#>        [1]     chr1       11874-12227      + |
#>        [2]     chr1       12613-12721      + |
#>        [3]     chr1       13221-14409      + |
#>        [4]     chr1       14362-14829      - |
#>        [5]     chr1       14970-15038      - |
#>        ...      ...               ...    ... .
#>   [459748]     chrY 59338754-59338859      + |
#>   [459749]     chrY 59338754-59338859      + |
#>   [459750]     chrY 59340194-59340278      + |
#>   [459751]     chrY 59342487-59343488      + |
#>   [459752]     chrY 59342487-59343488      + |
#>                                          name     score
#>                                   <character> <numeric>
#>        [1]    NR_046018_exon_0_0_chr1_11874_f         0
#>        [2]    NR_046018_exon_1_0_chr1_12613_f         0
#>        [3]    NR_046018_exon_2_0_chr1_13221_f         0
#>        [4]    NR_024540_exon_0_0_chr1_14362_r         0
#>        [5]    NR_024540_exon_1_0_chr1_14970_r         0
#>        ...                                ...       ...
#>   [459748] NM_002186_exon_6_0_chrY_59338754_f         0
#>   [459749] NM_176786_exon_7_0_chrY_59338754_f         0
#>   [459750] NM_002186_exon_7_0_chrY_59340194_f         0
#>   [459751] NM_002186_exon_8_0_chrY_59342487_f         0
#>   [459752] NM_176786_exon_8_0_chrY_59342487_f         0
#>   -------
#>   seqinfo: 93 sequences from an unspecified genome; no seqlengths
exons %>%
  group_by(seqnames) %>%
  summarise(n = n())
#> DataFrame with 49 rows and 2 columns
#>                 seqnames         n
#>                    <Rle> <integer>
#> 1                   chr1     43366
#> 2   chr1_gl000191_random        42
#> 3   chr1_gl000192_random        46
#> 4                  chr10     19420
#> 5                  chr11     24476
#> ...                  ...       ...
#> 45        chrUn_gl000222        20
#> 46        chrUn_gl000223        22
#> 47        chrUn_gl000228        85
#> 48                  chrX     18173
#> 49                  chrY      4128
```

Next we create a column representing the transcript\_id with `mutate`:

``` r
exons <- exons %>%
  mutate(tx_id = sub("_exon.*", "", name))
```

To find all GWAS SNPs that overlap exons, we use `join_overlap_inner`.
This will create a new *GRanges* with the coordinates of SNPs that
overlap exons, as well as metadata from both objects.

``` r
olap <- join_overlap_inner(gwas, exons)
olap
#> GRanges object with 3439 ranges and 4 metadata columns:
#>          seqnames    ranges strand |      name.x
#>             <Rle> <IRanges>  <Rle> | <character>
#>      [1]     chr1   1079198      * |  rs11260603
#>      [2]     chr1   1247494      * |     rs12103
#>      [3]     chr1   1247494      * |     rs12103
#>      [4]     chr1   1247494      * |     rs12103
#>      [5]     chr1   1247494      * |     rs12103
#>      ...      ...       ...    ... .         ...
#>   [3435]     chrX 153764217      * |   rs1050828
#>   [3436]     chrX 153764217      * |   rs1050828
#>   [3437]     chrX 153764217      * |   rs1050828
#>   [3438]     chrX 153764217      * |   rs1050828
#>   [3439]     chrX 153764217      * |   rs1050828
#>                                          name.y     score        tx_id
#>                                     <character> <numeric>  <character>
#>      [1]      NR_038869_exon_2_0_chr1_1078119_f         0    NR_038869
#>      [2]   NM_001256456_exon_1_0_chr1_1247398_r         0 NM_001256456
#>      [3]   NM_001256460_exon_1_0_chr1_1247398_r         0 NM_001256460
#>      [4]   NM_001256462_exon_1_0_chr1_1247398_r         0 NM_001256462
#>      [5]   NM_001256463_exon_1_0_chr1_1247398_r         0 NM_001256463
#>      ...                                    ...       ...          ...
#>   [3435] NM_001042351_exon_9_0_chrX_153764152_r         0 NM_001042351
#>   [3436]    NM_000402_exon_9_0_chrX_153764152_r         0    NM_000402
#>   [3437] NM_001042351_exon_9_0_chrX_153764152_r         0 NM_001042351
#>   [3438]    NM_000402_exon_9_0_chrX_153764152_r         0    NM_000402
#>   [3439] NM_001042351_exon_9_0_chrX_153764152_r         0 NM_001042351
#>   -------
#>   seqinfo: 93 sequences from an unspecified genome; no seqlengths
```

For each SNP we can count the number of times it overlaps a transcript.

``` r
olap %>%
  group_by(name.x, tx_id) %>%
  summarise(n = n())
#> DataFrame with 1619 rows and 3 columns
#>           name.x       tx_id         n
#>      <character> <character> <integer>
#> 1     rs17121403   NM_000028         1
#> 2       rs429358   NM_000041         3
#> 3         rs7412   NM_000041         3
#> 4      rs2234978   NM_000043         1
#> 5      rs1801516   NM_000051         1
#> ...          ...         ...       ...
#> 1615   rs1420101   NR_104167         1
#> 1616   rs1065656   NR_104318         1
#> 1617   rs2224391   NR_104417         1
#> 1618   rs2224391   NR_104418         1
#> 1619  rs41281112   NR_104592         1
```

We can also generate 2bp splice sites on either side of the exon using
`flank_left` and `flank_right`. We add a column indicating the side of
flanking for illustrative purposes. The `interweave` function pairs the
left and right ranges objects.

``` r
left_ss <- flank_left(exons, 2L) 
right_ss <- flank_right(exons, 2L)
all_ss <- interweave(left_ss, right_ss, .id = "side")
all_ss
#> GRanges object with 919504 ranges and 4 metadata columns:
#>            seqnames            ranges strand |
#>               <Rle>         <IRanges>  <Rle> |
#>        [1]     chr1       11872-11873      + |
#>        [2]     chr1       12228-12229      + |
#>        [3]     chr1       12611-12612      + |
#>        [4]     chr1       12722-12723      + |
#>        [5]     chr1       13219-13220      + |
#>        ...      ...               ...    ... .
#>   [919500]     chrY 59340279-59340280      + |
#>   [919501]     chrY 59342485-59342486      + |
#>   [919502]     chrY 59343489-59343490      + |
#>   [919503]     chrY 59342485-59342486      + |
#>   [919504]     chrY 59343489-59343490      + |
#>                                          name     score       tx_id
#>                                   <character> <numeric> <character>
#>        [1]    NR_046018_exon_0_0_chr1_11874_f         0   NR_046018
#>        [2]    NR_046018_exon_0_0_chr1_11874_f         0   NR_046018
#>        [3]    NR_046018_exon_1_0_chr1_12613_f         0   NR_046018
#>        [4]    NR_046018_exon_1_0_chr1_12613_f         0   NR_046018
#>        [5]    NR_046018_exon_2_0_chr1_13221_f         0   NR_046018
#>        ...                                ...       ...         ...
#>   [919500] NM_002186_exon_7_0_chrY_59340194_f         0   NM_002186
#>   [919501] NM_002186_exon_8_0_chrY_59342487_f         0   NM_002186
#>   [919502] NM_002186_exon_8_0_chrY_59342487_f         0   NM_002186
#>   [919503] NM_176786_exon_8_0_chrY_59342487_f         0   NM_176786
#>   [919504] NM_176786_exon_8_0_chrY_59342487_f         0   NM_176786
#>                   side
#>            <character>
#>        [1]        left
#>        [2]       right
#>        [3]        left
#>        [4]       right
#>        [5]        left
#>        ...         ...
#>   [919500]       right
#>   [919501]        left
#>   [919502]       right
#>   [919503]        left
#>   [919504]       right
#>   -------
#>   seqinfo: 93 sequences from an unspecified genome; no seqlengths
```

# Learning more

Read the [introduction
vignette](https://sa-lee.github.io/plyranges/articles/an-introduction.html)
for an overview of `plyranges` features.

# About the design of plyranges API

The `plyranges` package aims to provide a
[fluent](https://en.wikipedia.org/wiki/Fluent_interface) interface for
performing common genomic data analysis tasks. One goal of this package
is to decrease the learning curve for Bioconductor and S4 classes
(especially for new users) by providing a consistent interface to common
classes *IRanges* and *GRanges* (future work will extend this to the
*SummarizedExperiment* class)

All of the methods defined in `plyranges` are type-consistent and
chainable. As a consequence, users familiar with `dplyr` should be able
to read and understand `plyranges` code without difficulty.
