# Join data by metadata columns

These functions enable joins based on metadata of `x` and data in `y`,
utilizing
[`S4Vectors::merge()`](https://rdrr.io/pkg/S4Vectors/man/Vector-merge.html).

## Usage

``` r
join_mcols_inner(x, y, by = NULL, ...)

join_mcols_left(x, y, by = NULL, ...)
```

## Arguments

- x:

  Object representing ranges, with metadata columns containing variables
  for matching

- y:

  A table of data, DataFrame, data.frame or tibble

- by:

  Specifications of the columns used for merging. Passed to
  [`S4Vectors::merge()`](https://rdrr.io/pkg/S4Vectors/man/Vector-merge.html)

- ...:

  Additional arguments passed to
  [`S4Vectors::merge()`](https://rdrr.io/pkg/S4Vectors/man/Vector-merge.html)

## Value

Object representing ranges, with new metadata columns

## Details

The function `join_mcols_inner()` returns ranges in `x` only when a
match is found in the columns specified with `by` in the table `y`.

The function `join_mcols_left()` returns all ranges in `x` regardless of
a match in `y`, with duplications possible from multiple matches.

## Examples

``` r

x <- as_granges(data.frame(
  seqnames = "chr1", start = 1:5, width=10,
  id=c(1:4,4), foo=letters[1:5]
))
y <- DataFrame(id=c(1:2,2,4), bar=letters[6:9])

# metadata joins:
join_mcols_inner(x, y, by="id")
#> GRanges object with 5 ranges and 3 metadata columns:
#>       seqnames    ranges strand |        id         foo         bar
#>          <Rle> <IRanges>  <Rle> | <numeric> <character> <character>
#>   [1]     chr1      1-10      * |         1           a           f
#>   [2]     chr1      2-11      * |         2           b           g
#>   [3]     chr1      2-11      * |         2           b           h
#>   [4]     chr1      4-13      * |         4           d           i
#>   [5]     chr1      5-14      * |         4           e           i
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
join_mcols_left(x, y, by="id")
#> GRanges object with 6 ranges and 3 metadata columns:
#>       seqnames    ranges strand |        id         foo         bar
#>          <Rle> <IRanges>  <Rle> | <numeric> <character> <character>
#>   [1]     chr1      1-10      * |         1           a           f
#>   [2]     chr1      2-11      * |         2           b           g
#>   [3]     chr1      2-11      * |         2           b           h
#>   [4]     chr1      3-12      * |         3           c        <NA>
#>   [5]     chr1      4-13      * |         4           d           i
#>   [6]     chr1      5-14      * |         4           e           i
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
