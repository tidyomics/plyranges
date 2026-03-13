# Choose rows by their position

Choose rows by their position

## Usage

``` r
# S3 method for class 'Ranges'
slice(.data, ..., .preserve = FALSE)

# S3 method for class 'GroupedGenomicRanges'
slice(.data, ..., .preserve = FALSE)

# S3 method for class 'GroupedIntegerRanges'
slice(.data, ..., .preserve = FALSE)
```

## Arguments

- .data:

  a `Ranges` object

- ...:

  Integer row values indicating rows to keep. If `.data` has been
  grouped via
  [`group_by.GenomicRanges()`](https://tidyomics.github.io/plyranges/reference/group_by-ranges.md),
  then the positions are selected within each group.

- .preserve:

  when FALSE (the default) the grouping structure is recomputed,
  otherwise it is kept as is. Currently ignored.

## Value

a GRanges object

## Examples

``` r
df <- data.frame(start = 1:10,
                 width = 5,
                 seqnames = "seq1",
                 strand = sample(c("+", "-", "*"), 10, replace = TRUE),
                 gc = runif(10))
rng <- as_granges(df)
dplyr::slice(rng, 1:2)
#> GRanges object with 2 ranges and 1 metadata column:
#>       seqnames    ranges strand |        gc
#>          <Rle> <IRanges>  <Rle> | <numeric>
#>   [1]     seq1       1-5      - |  0.792625
#>   [2]     seq1       2-6      - |  0.614470
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
dplyr::slice(rng, -n())
#> GRanges object with 9 ranges and 1 metadata column:
#>       seqnames    ranges strand |        gc
#>          <Rle> <IRanges>  <Rle> | <numeric>
#>   [1]     seq1       1-5      - |  0.792625
#>   [2]     seq1       2-6      - |  0.614470
#>   [3]     seq1       3-7      - |  0.700621
#>   [4]     seq1       4-8      - |  0.107933
#>   [5]     seq1       5-9      * |  0.333866
#>   [6]     seq1      6-10      + |  0.324132
#>   [7]     seq1      7-11      * |  0.565426
#>   [8]     seq1      8-12      + |  0.203178
#>   [9]     seq1      9-13      + |  0.219705
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
dplyr::slice(rng, -5:-n())
#> GRanges object with 4 ranges and 1 metadata column:
#>       seqnames    ranges strand |        gc
#>          <Rle> <IRanges>  <Rle> | <numeric>
#>   [1]     seq1       1-5      - |  0.792625
#>   [2]     seq1       2-6      - |  0.614470
#>   [3]     seq1       3-7      - |  0.700621
#>   [4]     seq1       4-8      - |  0.107933
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

by_strand <- group_by(rng, strand)

# slice with group by finds positions within each group
dplyr::slice(by_strand, n())
#> GRanges object with 3 ranges and 1 metadata column:
#> Groups: strand [3]
#>       seqnames    ranges strand |        gc
#>          <Rle> <IRanges>  <Rle> | <numeric>
#>   [1]     seq1      7-11      * |  0.565426
#>   [2]     seq1      9-13      + |  0.219705
#>   [3]     seq1     10-14      - |  0.633115
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
dplyr::slice(by_strand, which.max(gc))
#> GRanges object with 3 ranges and 1 metadata column:
#> Groups: strand [3]
#>       seqnames    ranges strand |        gc
#>          <Rle> <IRanges>  <Rle> | <numeric>
#>   [1]     seq1       1-5      - |  0.792625
#>   [2]     seq1      6-10      + |  0.324132
#>   [3]     seq1      7-11      * |  0.565426
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

# if the index is beyond the number of groups slice are ignored
dplyr::slice(by_strand, 1:3)
#> GRanges object with 8 ranges and 1 metadata column:
#> Groups: strand [3]
#>       seqnames    ranges strand |        gc
#>          <Rle> <IRanges>  <Rle> | <numeric>
#>   [1]     seq1       1-5      - |  0.792625
#>   [2]     seq1       2-6      - |  0.614470
#>   [3]     seq1       3-7      - |  0.700621
#>   [4]     seq1       5-9      * |  0.333866
#>   [5]     seq1      6-10      + |  0.324132
#>   [6]     seq1      7-11      * |  0.565426
#>   [7]     seq1      8-12      + |  0.203178
#>   [8]     seq1      9-13      + |  0.219705
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

```
