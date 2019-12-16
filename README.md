
<!-- README.md is generated from README.Rmd. Please edit that file -->

# plyranges: fluent genomic data analysis <img id="plyranges_logo" src="man/figures/logo.png" align="right" width = "125" />

<!-- badges: start -->

[![Build
Status](https://travis-ci.org/sa-lee/plyranges.svg?branch=master)](https://travis-ci.org/sa-lee/plyranges)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/sa-lee/plyranges?branch=master&svg=true)](https://ci.appveyor.com/project/sa-lee/plyranges)
[![Codecov test
coverage](https://codecov.io/gh/sa-lee/plyranges/branch/master/graph/badge.svg)](https://codecov.io/gh/sa-lee/plyranges?branch=master)
<!-- badges: end -->

[plyranges](https://www.bioconductor.org/packages/release/bioc/html/plyranges.html)
provides a consistent interface for importing and wrangling genomics
data from a variety of sources. The package defines a grammar of genomic
data transformation based on `dplyr` and the Bioconductor packages
`IRanges`, `GenomicRanges`, and `rtracklayer`. It does this by providing
a set of verbs for developing analysis pipelines based on *Ranges*
objects that represent genomic regions:

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

For more details on the features of plyranges, read the
[vignette](https://sa-lee.github.io/plyranges/articles/an-introduction.html).
For a complete case-study on using plyranges to combine ATAC-seq and
RNA-seq results read the [*fluentGenomics*
workflow](https://sa-lee.github.io/fluentGenomics).

# Installation

[plyranges](https://www.bioconductor.org/packages/release/bioc/html/plyranges.html)
can be installed from the latest Bioconductor release:

``` r
# install.packages("BiocManager")
BiocManager::install("plyranges")
```

To install the development version from GitHub:

``` r
BiocManager::install("sa-lee/plyranges")
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
#>           name.x        tx_id         n
#>      <character>  <character> <integer>
#> 1     rs10043775 NM_001271723         1
#> 2     rs10043775    NM_030793         1
#> 3        rs10078 NM_001242412         1
#> 4        rs10078    NM_020731         1
#> 5        rs10089    NM_001046         1
#> ...          ...          ...       ...
#> 1615   rs9906595 NM_001008777         1
#> 1616      rs9948    NM_017623         1
#> 1617      rs9948    NM_199078         1
#> 1618    rs995030    NM_000899         4
#> 1619    rs995030    NM_003994         4
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
#>                                          name     score       tx_id        side
#>                                   <character> <numeric> <character> <character>
#>        [1]    NR_046018_exon_0_0_chr1_11874_f         0   NR_046018        left
#>        [2]    NR_046018_exon_0_0_chr1_11874_f         0   NR_046018       right
#>        [3]    NR_046018_exon_1_0_chr1_12613_f         0   NR_046018        left
#>        [4]    NR_046018_exon_1_0_chr1_12613_f         0   NR_046018       right
#>        [5]    NR_046018_exon_2_0_chr1_13221_f         0   NR_046018        left
#>        ...                                ...       ...         ...         ...
#>   [919500] NM_002186_exon_7_0_chrY_59340194_f         0   NM_002186       right
#>   [919501] NM_002186_exon_8_0_chrY_59342487_f         0   NM_002186        left
#>   [919502] NM_002186_exon_8_0_chrY_59342487_f         0   NM_002186       right
#>   [919503] NM_176786_exon_8_0_chrY_59342487_f         0   NM_176786        left
#>   [919504] NM_176786_exon_8_0_chrY_59342487_f         0   NM_176786       right
#>   -------
#>   seqinfo: 93 sequences from an unspecified genome; no seqlengths
```

# Learning more

  - The [*fluentGenomics*
    workflow](https://sa-lee.github.io/fluentGenomics) package shows you
    how to combine differential expression genes and differential
    chromatin accessibility peaks using plyranges. It extends the [case
    study](https://github.com/mikelove/plyrangesTximetaCaseStudy) by
    Michael Love for using plyranges with
    [tximeta](https://bioconductor.org/packages/release/bioc/html/tximeta.html).

  - The [extended vignette in the plyrangesWorkshops
    package](https://github.com/sa-lee/plyrangesWorkshops) has a
    detailed walk through of using plyranges for coverage analysis.

  - The [Bioc 2018 Workshop
    book](https://bioconductor.github.io/BiocWorkshops/fluent-genomic-data-analysis-with-plyranges.html)
    has worked examples of using `plyranges` to analyse publicly
    available genomics data.

# Citation

If you found `plyranges` useful for your work please cite our
[paper](http://dx.doi.org/10.1186/s13059-018-1597-8):

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

# Contributing

We welcome contributions from the R/Bioconductor community. We ask that
contributors follow the [code of conduct](CODE_OF_CONDUCT.md) and the
guide outlined [here](CONTRIBUTING.md).
