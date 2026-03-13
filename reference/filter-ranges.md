# Subset a `Ranges` object

Subset a `Ranges` object

## Usage

``` r
# S3 method for class 'Ranges'
filter(.data, ..., .preserve = FALSE)
```

## Arguments

- .data:

  A `Ranges` object

- ...:

  valid logical predictates to subset .data by. These are determined by
  variables in `.data`. If more than one condition is supplied, the
  conditions are combined with `&`. Only rows where the condition
  evaluates to `TRUE` are kept.

- .preserve:

  when FALSE (the default) grouping structure is recalculated, TRUE is
  currently not implemented.

## Value

a Ranges object

## Details

For any Ranges objects `filter` can act on all core components of the
class including start, end, width (for IRanges) or seqnames and strand
(for GRanges) in addition to metadata columns. If the Ranges object is
grouped, `filter` will act seperately on each parition of the data.

## See also

[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)

## Examples

``` r
set.seed(100)
df <- data.frame(start = 1:10,
                 width = 5,
                 seqnames = "seq1",
                 strand = sample(c("+", "-", "*"), 10, replace = TRUE),
                 gc = runif(10))

rng <- as_granges(df)

filter(rng, strand == "+")
#> GRanges object with 1 range and 1 metadata column:
#>       seqnames    ranges strand |        gc
#>          <Rle> <IRanges>  <Rle> | <numeric>
#>   [1]     seq1       5-9      + |  0.357525
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
filter(rng, gc > 0.5)
#> GRanges object with 6 ranges and 1 metadata column:
#>       seqnames    ranges strand |        gc
#>          <Rle> <IRanges>  <Rle> | <numeric>
#>   [1]     seq1       2-6      * |  0.762551
#>   [2]     seq1       3-7      - |  0.669022
#>   [3]     seq1      7-11      - |  0.690291
#>   [4]     seq1      8-12      * |  0.535811
#>   [5]     seq1      9-13      - |  0.710804
#>   [6]     seq1     10-14      - |  0.538349
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

# multiple criteria
filter(rng, strand == "+" | start > 5)
#> GRanges object with 6 ranges and 1 metadata column:
#>       seqnames    ranges strand |        gc
#>          <Rle> <IRanges>  <Rle> | <numeric>
#>   [1]     seq1       5-9      + |  0.357525
#>   [2]     seq1      6-10      - |  0.359475
#>   [3]     seq1      7-11      - |  0.690291
#>   [4]     seq1      8-12      * |  0.535811
#>   [5]     seq1      9-13      - |  0.710804
#>   [6]     seq1     10-14      - |  0.538349
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
filter(rng, strand == "+" & start > 5)
#> GRanges object with 0 ranges and 1 metadata column:
#>    seqnames    ranges strand |        gc
#>       <Rle> <IRanges>  <Rle> | <numeric>
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

# multiple conditions are the same as and
filter(rng, strand == "+", start > 5)
#> GRanges object with 0 ranges and 1 metadata column:
#>    seqnames    ranges strand |        gc
#>       <Rle> <IRanges>  <Rle> | <numeric>
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

# grouping acts on each subset of the data
rng %>%
  group_by(strand) %>%
  filter(gc > 0.5)
#> GRanges object with 6 ranges and 1 metadata column:
#> Groups: strand [2]
#>       seqnames    ranges strand |        gc
#>          <Rle> <IRanges>  <Rle> | <numeric>
#>   [1]     seq1       2-6      * |  0.762551
#>   [2]     seq1       3-7      - |  0.669022
#>   [3]     seq1      7-11      - |  0.690291
#>   [4]     seq1      8-12      * |  0.535811
#>   [5]     seq1      9-13      - |  0.710804
#>   [6]     seq1     10-14      - |  0.538349
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
