# Shift all coordinates in a genomic interval left or right, upstream or downstream

Shift all coordinates in a genomic interval left or right, upstream or
downstream

## Usage

``` r
shift_left(x, shift = 0L)

shift_right(x, shift = 0L)

shift_upstream(x, shift = 0L)

shift_downstream(x, shift = 0L)
```

## Arguments

- x:

  a Ranges object .

- shift:

  the amount to move the genomic interval in the Ranges object by.
  Either a non-negative integer vector of length 1 or an integer vector
  the same length as x.

## Value

a Ranges object with start and end coordinates shifted.

## Details

Shifting left or right will ignore any strand information in the Ranges
object, while shifting upstream/downstream will shift coordinates on the
positive strand left/right and the negative strand right/left. By
default, unstranded features are treated as positive. When using
`shift_upstream()` or `shift_downstream()` when the `shift` argument is
indexed by the strandedness of the input ranges.

## See also

`IRanges::`[`shift()`](https://rdrr.io/pkg/IRanges/man/intra-range-methods.html),
`GenomicRanges::`[`shift()`](https://rdrr.io/pkg/GenomicRanges/man/intra-range-methods.html)

## Examples

``` r
ir <- as_iranges(data.frame(start = 10:15, width = 5))
shift_left(ir, 5L)
#> IRanges object with 6 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]         5         9         5
#>   [2]         6        10         5
#>   [3]         7        11         5
#>   [4]         8        12         5
#>   [5]         9        13         5
#>   [6]        10        14         5
shift_right(ir, 5L)
#> IRanges object with 6 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]        15        19         5
#>   [2]        16        20         5
#>   [3]        17        21         5
#>   [4]        18        22         5
#>   [5]        19        23         5
#>   [6]        20        24         5
gr <- as_granges(data.frame(start = 10:15,
                            width = 5,
                            seqnames = "seq1",
                            strand = c("+", "+", "-", "-", "+", "*")))
shift_upstream(gr, 5L)
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
shift_downstream(gr, 5L)
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
