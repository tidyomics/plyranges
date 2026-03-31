# Choose ranges by their position

The `slice` family of functions let you select ranges by position. They
follow the same conventions as their dplyr counterparts for data frames:

- `slice()` selects ranges by explicit integer position.

- `slice_head()` / `slice_tail()` select the first or last ranges.

- `slice_sample()` randomly selects ranges.

- `slice_min()` / `slice_max()` select ranges with the smallest or
  largest values of a metadata column.

All functions preserve grouping: when `.data` has been grouped with
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html),
selection is performed within each group.

## Usage

``` r
# S3 method for class 'Ranges'
slice(.data, ..., .preserve = FALSE)

# S3 method for class 'GroupedGenomicRanges'
slice(.data, ..., .preserve = FALSE)

# S3 method for class 'GroupedIntegerRanges'
slice(.data, ..., .preserve = FALSE)

# S3 method for class 'Ranges'
slice_head(.data, ..., n = NULL, prop = NULL, .preserve = FALSE)

# S3 method for class 'GroupedGenomicRanges'
slice_head(.data, ..., n = NULL, prop = NULL, .preserve = FALSE)

# S3 method for class 'GroupedIntegerRanges'
slice_head(.data, ..., n = NULL, prop = NULL, .preserve = FALSE)

# S3 method for class 'Ranges'
slice_tail(.data, ..., n = NULL, prop = NULL, .preserve = FALSE)

# S3 method for class 'GroupedGenomicRanges'
slice_tail(.data, ..., n = NULL, prop = NULL, .preserve = FALSE)

# S3 method for class 'GroupedIntegerRanges'
slice_tail(.data, ..., n = NULL, prop = NULL, .preserve = FALSE)

# S3 method for class 'Ranges'
slice_sample(
  .data,
  ...,
  n = NULL,
  prop = NULL,
  weight_by = NULL,
  replace = FALSE,
  .preserve = FALSE
)

# S3 method for class 'GroupedGenomicRanges'
slice_sample(
  .data,
  ...,
  n = NULL,
  prop = NULL,
  weight_by = NULL,
  replace = FALSE,
  .preserve = FALSE
)

# S3 method for class 'GroupedIntegerRanges'
slice_sample(
  .data,
  ...,
  n = NULL,
  prop = NULL,
  weight_by = NULL,
  replace = FALSE,
  .preserve = FALSE
)

# S3 method for class 'Ranges'
slice_min(
  .data,
  order_by,
  ...,
  n = NULL,
  prop = NULL,
  with_ties = TRUE,
  .preserve = FALSE
)

# S3 method for class 'GroupedGenomicRanges'
slice_min(
  .data,
  order_by,
  ...,
  n = NULL,
  prop = NULL,
  with_ties = TRUE,
  .preserve = FALSE
)

# S3 method for class 'GroupedIntegerRanges'
slice_min(
  .data,
  order_by,
  ...,
  n = NULL,
  prop = NULL,
  with_ties = TRUE,
  .preserve = FALSE
)

# S3 method for class 'Ranges'
slice_max(
  .data,
  order_by,
  ...,
  n = NULL,
  prop = NULL,
  with_ties = TRUE,
  .preserve = FALSE
)

# S3 method for class 'GroupedGenomicRanges'
slice_max(
  .data,
  order_by,
  ...,
  n = NULL,
  prop = NULL,
  with_ties = TRUE,
  .preserve = FALSE
)

# S3 method for class 'GroupedIntegerRanges'
slice_max(
  .data,
  order_by,
  ...,
  n = NULL,
  prop = NULL,
  with_ties = TRUE,
  .preserve = FALSE
)
```

## Arguments

- .data:

  a `Ranges` object.

- ...:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Integer row values indicating positions to keep. Only used by
  `slice()`; ignored elsewhere.

- .preserve:

  Ignored; retained for compatibility with the dplyr generic.

- n:

  Number of ranges to select. Mutually exclusive with `prop`. For
  `slice_min()` and `slice_max()`, negative values are not supported.

- prop:

  Fraction of ranges to select, in `(0, 1]`. Mutually exclusive with
  `n`.

- weight_by:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Unquoted name of a non-negative numeric metadata column used as
  sampling weights. Used by `slice_sample()`.

- replace:

  Should sampling be with replacement? Default `FALSE`. Used by
  `slice_sample()`.

- order_by:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Unquoted name of a numeric metadata column to order by. Used by
  `slice_min()` and `slice_max()`.

- with_ties:

  If `TRUE` (the default), all ranges tied at the boundary rank are
  kept, which may return more than `n` ranges. Set to `FALSE` to return
  exactly `n`. Used by `slice_min()` and `slice_max()`.

## Value

A `GRanges` (or `IRanges`) object of the same class as `.data`, with
grouping structure re-applied when `.data` was grouped.

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
