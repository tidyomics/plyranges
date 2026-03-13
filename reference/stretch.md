# Stretch a genomic interval

By default, `stretch(x)` will anchor by the center of a Ranges object.
This means that half of the value of `extend` will be added to the end
of the range and the remaining half subtracted from the start of the
Range. The other anchors will leave the start/end fixed and stretch the
end/start respectively.

## Usage

``` r
stretch(x, extend)
```

## Arguments

- x:

  a Ranges object, to fix by either the start, end or center of an
  interval use `anchor_start(x)`, `anchor_end(x)`, `anchor_center(x)`.
  To fix by strand use `anchor_3p(x)` or `anchor_5p(x)`.

- extend:

  the amount to alter the width of a Ranges object by. Either an integer
  vector of length 1 or an integer vector the same length as x.

## Value

a Ranges object with modified start or end (or both) coordinates

## See also

[`anchor()`](https://tidyomics.github.io/plyranges/reference/ranges-anchor.md),
[`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html)

## Examples

``` r
rng <- as_iranges(data.frame(start=c(2:-1, 13:15), width=c(0:3, 2:0)))
rng2 <- stretch(anchor_center(rng), 10)
stretch(anchor_start(rng2), 10)
#> IRanges object with 7 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]        -3        16        20
#>   [2]        -4        16        21
#>   [3]        -5        16        22
#>   [4]        -6        16        23
#>   [5]         8        29        22
#>   [6]         9        29        21
#>   [7]        10        29        20
stretch(anchor_end(rng2), 10)
#> IRanges object with 7 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]       -13         6        20
#>   [2]       -14         6        21
#>   [3]       -15         6        22
#>   [4]       -16         6        23
#>   [5]        -2        19        22
#>   [6]        -1        19        21
#>   [7]         0        19        20
grng <- as_granges(data.frame(seqnames = "chr1",
                         strand = c("+", "-", "-", "+", "+", "-", "+"),
                         start=c(2:-1, 13:15),
                         width=c(0:3, 2:0)))
stretch(anchor_3p(grng), 10)
#> GRanges object with 7 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1      -8-1      +
#>   [2]     chr1      1-11      -
#>   [3]     chr1      0-11      -
#>   [4]     chr1     -11-1      +
#>   [5]     chr1      3-14      +
#>   [6]     chr1     14-24      -
#>   [7]     chr1      5-14      +
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
stretch(anchor_5p(grng), 10)
#> GRanges object with 7 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1      2-11      +
#>   [2]     chr1      -9-1      -
#>   [3]     chr1     -10-1      -
#>   [4]     chr1     -1-11      +
#>   [5]     chr1     13-24      +
#>   [6]     chr1      4-14      -
#>   [7]     chr1     15-24      +
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
