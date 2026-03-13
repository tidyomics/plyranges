# Generate flanking regions

Find flanking regions to the left or right or upstream or downstream of
a Ranges object.

## Usage

``` r
flank_left(x, width = 0L)

flank_right(x, width = 0L)

flank_upstream(x, width = 0L)

flank_downstream(x, width = 0L)
```

## Arguments

- x:

  a Ranges object.

- width:

  the width of the flanking region relative to the ranges in `x`. Either
  an integer vector of length 1 or an integer vector the same length
  as x. The width can be negative in which case the flanking region is
  reversed.

## Value

A Ranges object of same length as `x`.

## Details

The function `flank_left` will create the flanking region to the left of
starting coordinates in `x`, while `flank_right` will create the
flanking region to the right of the starting coordinates in `x`. The
function `flank_upstream` will `flank_left` if the strand of rows in `x`
is not negative and will `flank_right` if the strand of rows in `x` is
negative. The function `flank_downstream` will `flank_right` if the
strand of rows in `x` is not negative and will `flank_leftt` if the
strand of rows in `x` is negative.

By default `flank_left` and `flank_right` will ignore strandedness of
any ranges, while `flank_upstream` and `flank_downstream` will take into
account the strand of `x`.

## See also

`IRanges::`[`flank()`](https://rdrr.io/pkg/IRanges/man/intra-range-methods.html),
`GenomicRanges::`[`flank()`](https://rdrr.io/pkg/GenomicRanges/man/intra-range-methods.html)

## Examples

``` r
gr <- as_granges(data.frame(start = 10:15,
                            width = 5,
                            seqnames = "seq1",
                            strand = c("+", "+", "-", "-", "+", "*")))
flank_left(gr, width = 5L)
#> GRanges object with 6 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     seq1       5-9      +
#>   [2]     seq1      6-10      +
#>   [3]     seq1      7-11      -
#>   [4]     seq1      8-12      -
#>   [5]     seq1      9-13      +
#>   [6]     seq1     10-14      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
flank_right(gr, width = 5L)
#> GRanges object with 6 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     seq1     15-19      +
#>   [2]     seq1     16-20      +
#>   [3]     seq1     17-21      -
#>   [4]     seq1     18-22      -
#>   [5]     seq1     19-23      +
#>   [6]     seq1     20-24      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
flank_upstream(gr, width = 5L)
#> GRanges object with 6 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     seq1       5-9      +
#>   [2]     seq1      6-10      +
#>   [3]     seq1     17-21      -
#>   [4]     seq1     18-22      -
#>   [5]     seq1      9-13      +
#>   [6]     seq1     10-14      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
flank_downstream(gr, width = 5L)
#> GRanges object with 6 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     seq1     15-19      +
#>   [2]     seq1     16-20      +
#>   [3]     seq1      7-11      -
#>   [4]     seq1      8-12      -
#>   [5]     seq1     19-23      +
#>   [6]     seq1     20-24      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
