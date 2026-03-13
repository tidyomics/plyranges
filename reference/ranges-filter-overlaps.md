# Filter by overlapping/non-overlapping ranges

Filter by overlapping/non-overlapping ranges

## Usage

``` r
filter_by_overlaps(x, y, maxgap = -1L, minoverlap = 0L)

filter_by_non_overlaps(x, y, maxgap, minoverlap)

filter_by_overlaps_directed(x, y, maxgap = -1L, minoverlap = 0L)

filter_by_non_overlaps_directed(x, y, maxgap, minoverlap)
```

## Arguments

- x, y:

  Objects representing ranges

- maxgap:

  The maximimum gap between intervals as a single integer greater than
  or equal to -1. If you modify this argument, `minoverlap` must be held
  fixed.

- minoverlap:

  The minimum amount of overlap between intervals as a single integer
  greater than 0. If you modify this argument, `maxgap` must be held
  fixed.

## Value

a Ranges object

## Details

By default, `filter_by_overlaps` and `filter_by_non_overlaps` ignore
strandedness for
[`GenomicRanges::GRanges()`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
objects. To perform stranded operations use
`filter_by_overlaps_directed` and `filter_by_non_overlaps_directed`. The
argument `maxgap` is the maximum number of positions between two ranges
for them to be considered overlapping. Here the default is set to be -1
as that is the the gap between two ranges that has its start or end
strictly inside the other. The argugment `minoverlap` refers to the
minimum number of positions overlapping between ranges, to consider
there to be overlap.

## See also

`IRanges::`[`subsetByOverlaps()`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html)

## Examples

``` r
df <- data.frame(seqnames = c("chr1", rep("chr2", 2),
                              rep("chr3", 3), rep("chr4", 4)),
                 start = 1:10,
                 width = 10:1,
                 strand = c("-", "+", "+", "*", "*", "+", "+", "+", "-", "-"),
                 name = letters[1:10])
query <- as_granges(df)

df2 <- data.frame(seqnames = c(rep("chr2", 2), rep("chr1", 3), "chr2"),
                  start = c(4,3,7,13,1,4),
                  width = c(6,6,3,3,3,9),
                  strand = c(rep("+", 3), rep("-", 3)))
subject <- as_granges(df2)

filter_by_overlaps(query, subject)
#> GRanges object with 3 ranges and 1 metadata column:
#>       seqnames    ranges strand |        name
#>          <Rle> <IRanges>  <Rle> | <character>
#>   [1]     chr1      1-10      - |           a
#>   [2]     chr2      2-10      + |           b
#>   [3]     chr2      3-10      + |           c
#>   -------
#>   seqinfo: 4 sequences from an unspecified genome; no seqlengths

filter_by_overlaps_directed(query, subject)
#> GRanges object with 3 ranges and 1 metadata column:
#>       seqnames    ranges strand |        name
#>          <Rle> <IRanges>  <Rle> | <character>
#>   [1]     chr1      1-10      - |           a
#>   [2]     chr2      2-10      + |           b
#>   [3]     chr2      3-10      + |           c
#>   -------
#>   seqinfo: 4 sequences from an unspecified genome; no seqlengths

filter_by_non_overlaps(query, subject)
#> GRanges object with 7 ranges and 1 metadata column:
#>       seqnames    ranges strand |        name
#>          <Rle> <IRanges>  <Rle> | <character>
#>   [1]     chr3      4-10      * |           d
#>   [2]     chr3      5-10      * |           e
#>   [3]     chr3      6-10      + |           f
#>   [4]     chr4      7-10      + |           g
#>   [5]     chr4      8-10      + |           h
#>   [6]     chr4      9-10      - |           i
#>   [7]     chr4        10      - |           j
#>   -------
#>   seqinfo: 4 sequences from an unspecified genome; no seqlengths

filter_by_non_overlaps_directed(query, subject)
#> GRanges object with 7 ranges and 1 metadata column:
#>       seqnames    ranges strand |        name
#>          <Rle> <IRanges>  <Rle> | <character>
#>   [1]     chr3      4-10      * |           d
#>   [2]     chr3      5-10      * |           e
#>   [3]     chr3      6-10      + |           f
#>   [4]     chr4      7-10      + |           g
#>   [5]     chr4      8-10      + |           h
#>   [6]     chr4      9-10      - |           i
#>   [7]     chr4        10      - |           j
#>   -------
#>   seqinfo: 4 sequences from an unspecified genome; no seqlengths
```
