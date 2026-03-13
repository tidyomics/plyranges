# Interweave a pair of Ranges objects together

Interweave a pair of Ranges objects together

## Usage

``` r
interweave(left, right, .id = NULL)
```

## Arguments

- left, right:

  Ranges objects.

- .id:

  When supplied a new column that represents the origin column and is
  linked to each row of the resulting Ranges object.

## Value

a Ranges object

## Details

The output of `interweave()` takes pairs of Ranges objects and combines
them into a single Ranges object. If an .id argument is supplied, an
origin column with name .id is created indicated which side the
resulting Range comes from (eit)

## Examples

``` r
gr <- as_granges(data.frame(start = 10:15,
                            width = 5,
                            seqnames = "seq1",
                            strand = c("+", "+", "-", "-", "+", "*")))
interweave(flank_left(gr, width = 5L), flank_right(gr, width = 5L))
#> GRanges object with 12 ranges and 0 metadata columns:
#>        seqnames    ranges strand
#>           <Rle> <IRanges>  <Rle>
#>    [1]     seq1       5-9      +
#>    [2]     seq1     15-19      +
#>    [3]     seq1      6-10      +
#>    [4]     seq1     16-20      +
#>    [5]     seq1      7-11      -
#>    ...      ...       ...    ...
#>    [8]     seq1     18-22      -
#>    [9]     seq1      9-13      +
#>   [10]     seq1     19-23      +
#>   [11]     seq1     10-14      *
#>   [12]     seq1     20-24      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
interweave(flank_left(gr, width = 5L), flank_right(gr, width = 5L), .id = "origin")
#> GRanges object with 12 ranges and 1 metadata column:
#>        seqnames    ranges strand |      origin
#>           <Rle> <IRanges>  <Rle> | <character>
#>    [1]     seq1       5-9      + |        left
#>    [2]     seq1     15-19      + |       right
#>    [3]     seq1      6-10      + |        left
#>    [4]     seq1     16-20      + |       right
#>    [5]     seq1      7-11      - |        left
#>    ...      ...       ...    ... .         ...
#>    [8]     seq1     18-22      - |       right
#>    [9]     seq1      9-13      + |        left
#>   [10]     seq1     19-23      + |       right
#>   [11]     seq1     10-14      * |        left
#>   [12]     seq1     20-24      * |       right
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
