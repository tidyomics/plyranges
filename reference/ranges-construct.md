# Construct a I/GRanges object from a tibble or data.frame

The as_i(g)ranges function looks for column names in .data called start,
end, width, seqnames and strand in order to construct an IRanges or
GRanges object. By default other columns in .data are placed into the
mcols ( metadata columns) slot of the returned object.

## Usage

``` r
as_iranges(.data, ..., keep_mcols = TRUE)

as_granges(.data, ..., keep_mcols = TRUE)
```

## Arguments

- .data:

  a [`data.frame()`](https://rdrr.io/r/base/data.frame.html) or
  [`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)
  to construct a Ranges object from

- ...:

  optional named arguments specifying which the columns in .data
  containin the core components a Ranges object.

- keep_mcols:

  place the remaining columns into the metadata columns slot
  (default=TRUE)

## Value

a Ranges object.

## See also

`IRanges::`[`IRanges()`](https://rdrr.io/pkg/IRanges/man/IRanges-constructor.html),
`GenomicRanges::`[`GRanges()`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)

## Examples

``` r
df <- data.frame(start=c(2:-1, 13:15), width=c(0:3, 2:0))
as_iranges(df)
#> IRanges object with 7 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]         2         1         0
#>   [2]         1         1         1
#>   [3]         0         1         2
#>   [4]        -1         1         3
#>   [5]        13        14         2
#>   [6]        14        14         1
#>   [7]        15        14         0

df <- data.frame(start=c(2:-1, 13:15), width=c(0:3, 2:0), strand = "+")
# will return an IRanges object
as_iranges(df)
#> IRanges object with 7 ranges and 1 metadata column:
#>           start       end     width |      strand
#>       <integer> <integer> <integer> | <character>
#>   [1]         2         1         0 |           +
#>   [2]         1         1         1 |           +
#>   [3]         0         1         2 |           +
#>   [4]        -1         1         3 |           +
#>   [5]        13        14         2 |           +
#>   [6]        14        14         1 |           +
#>   [7]        15        14         0 |           +

df <- data.frame(start=c(2:-1, 13:15), width=c(0:3, 2:0),
strand = "+", seqnames = "chr1")
as_granges(df)
#> GRanges object with 7 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1       2-1      +
#>   [2]     chr1         1      +
#>   [3]     chr1       0-1      +
#>   [4]     chr1      -1-1      +
#>   [5]     chr1     13-14      +
#>   [6]     chr1        14      +
#>   [7]     chr1     15-14      +
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

# as_g/iranges understand alternate name specification
df <- data.frame(start=c(2:-1, 13:15), width=c(0:3, 2:0),
strand = "+", chr = "chr1")
as_granges(df, seqnames = chr)
#> GRanges object with 7 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1       2-1      +
#>   [2]     chr1         1      +
#>   [3]     chr1       0-1      +
#>   [4]     chr1      -1-1      +
#>   [5]     chr1     13-14      +
#>   [6]     chr1        14      +
#>   [7]     chr1     15-14      +
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

# can also handle DFrame input
df <- methods::as(df, "DFrame")
df$y <- IRanges::IntegerList(c(1,2,3), NA, 5, 6, 8, 9, 10:12)
as_iranges(df)
#> IRanges object with 7 ranges and 3 metadata columns:
#>           start       end     width |      strand         chr             y
#>       <integer> <integer> <integer> | <character> <character> <IntegerList>
#>   [1]         2         1         0 |           +        chr1         1,2,3
#>   [2]         1         1         1 |           +        chr1          <NA>
#>   [3]         0         1         2 |           +        chr1             5
#>   [4]        -1         1         3 |           +        chr1             6
#>   [5]        13        14         2 |           +        chr1             8
#>   [6]        14        14         1 |           +        chr1             9
#>   [7]        15        14         0 |           +        chr1      10,11,12
as_granges(df, seqnames = chr)
#> GRanges object with 7 ranges and 1 metadata column:
#>       seqnames    ranges strand |             y
#>          <Rle> <IRanges>  <Rle> | <IntegerList>
#>   [1]     chr1       2-1      + |         1,2,3
#>   [2]     chr1         1      + |          <NA>
#>   [3]     chr1       0-1      + |             5
#>   [4]     chr1      -1-1      + |             6
#>   [5]     chr1     13-14      + |             8
#>   [6]     chr1        14      + |             9
#>   [7]     chr1     15-14      + |      10,11,12
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
