# Combine Ranges by concatentating them together

Combine Ranges by concatentating them together

## Usage

``` r
bind_ranges(..., .id = NULL)
```

## Arguments

- ...:

  Ranges objects to combine. Each argument can be a Ranges object, or a
  list of Ranges objects.

- .id:

  Ranges object identifier. When .id is supplied a new column is created
  that links each row to the original Range object. The contents of the
  column correspond to the named arguments or the names of the list
  supplied.

## Value

a concatenated Ranges object

## Note

Currently GRangesList or IRangesList objects are not supported.

## Examples

``` r
gr <- as_granges(data.frame(start = 10:15,
                            width = 5,
                            seqnames = "seq1"))
gr2 <- as_granges(data.frame(start = 11:14,
                            width = 1:4,
                            seqnames = "seq2"))
bind_ranges(gr, gr2)
#> GRanges object with 10 ranges and 0 metadata columns:
#>        seqnames    ranges strand
#>           <Rle> <IRanges>  <Rle>
#>    [1]     seq1     10-14      *
#>    [2]     seq1     11-15      *
#>    [3]     seq1     12-16      *
#>    [4]     seq1     13-17      *
#>    [5]     seq1     14-18      *
#>    [6]     seq1     15-19      *
#>    [7]     seq2        11      *
#>    [8]     seq2     12-13      *
#>    [9]     seq2     13-15      *
#>   [10]     seq2     14-17      *
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

bind_ranges(a = gr, b = gr2, .id = "origin")
#> GRanges object with 10 ranges and 1 metadata column:
#>        seqnames    ranges strand | origin
#>           <Rle> <IRanges>  <Rle> |  <Rle>
#>    [1]     seq1     10-14      * |      a
#>    [2]     seq1     11-15      * |      a
#>    [3]     seq1     12-16      * |      a
#>    [4]     seq1     13-17      * |      a
#>    [5]     seq1     14-18      * |      a
#>    [6]     seq1     15-19      * |      a
#>    [7]     seq2        11      * |      b
#>    [8]     seq2     12-13      * |      b
#>    [9]     seq2     13-15      * |      b
#>   [10]     seq2     14-17      * |      b
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

bind_ranges(gr, list(gr, gr2), gr2)
#> GRanges object with 20 ranges and 0 metadata columns:
#>        seqnames    ranges strand
#>           <Rle> <IRanges>  <Rle>
#>    [1]     seq1     10-14      *
#>    [2]     seq1     11-15      *
#>    [3]     seq1     12-16      *
#>    [4]     seq1     13-17      *
#>    [5]     seq1     14-18      *
#>    ...      ...       ...    ...
#>   [16]     seq2     14-17      *
#>   [17]     seq2        11      *
#>   [18]     seq2     12-13      *
#>   [19]     seq2     13-15      *
#>   [20]     seq2     14-17      *
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

bind_ranges(list(a = gr, b = gr2), c = gr, .id = "origin")
#> GRanges object with 16 ranges and 1 metadata column:
#>        seqnames    ranges strand | origin
#>           <Rle> <IRanges>  <Rle> |  <Rle>
#>    [1]     seq1     10-14      * |      a
#>    [2]     seq1     11-15      * |      a
#>    [3]     seq1     12-16      * |      a
#>    [4]     seq1     13-17      * |      a
#>    [5]     seq1     14-18      * |      a
#>    ...      ...       ...    ... .    ...
#>   [12]     seq1     11-15      * |      c
#>   [13]     seq1     12-16      * |      c
#>   [14]     seq1     13-17      * |      c
#>   [15]     seq1     14-18      * |      c
#>   [16]     seq1     15-19      * |      c
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```
