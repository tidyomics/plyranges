# Coerce an Rle or RleList object to Ranges

Coerce an Rle or RleList object to Ranges

## Usage

``` r
as_ranges(.data)
```

## Arguments

- .data:

  a
  [`S4Vectors::Rle()`](https://rdrr.io/pkg/S4Vectors/man/Rle-class.html)
  or an
  [`IRanges::RleList()`](https://rdrr.io/pkg/IRanges/man/AtomicList-class.html)
  object.

## Value

an
[`IRanges::IRanges()`](https://rdrr.io/pkg/IRanges/man/IRanges-constructor.html)
object if the input is an
[`S4Vectors::Rle()`](https://rdrr.io/pkg/S4Vectors/man/Rle-class.html)
object or a
[`GenomicRanges::GRanges()`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object for an
[`IRanges::RleList()`](https://rdrr.io/pkg/IRanges/man/AtomicList-class.html)
object.

## Details

This function is behind
[`compute_coverage()`](https://tidyomics.github.io/plyranges/reference/compute_coverage.md).

## See also

`S4Vectors::`[`Rle()`](https://rdrr.io/pkg/S4Vectors/man/Rle-class.html),
`IRanges::`[`RleList()`](https://rdrr.io/pkg/IRanges/man/AtomicList-class.html)

## Examples

``` r
x <- S4Vectors::Rle(10:1, 1:10)
as_ranges(x)
#> IRanges object with 10 ranges and 1 metadata column:
#>            start       end     width |     score
#>        <integer> <integer> <integer> | <integer>
#>    [1]         1         1         1 |        10
#>    [2]         2         3         2 |         9
#>    [3]         4         6         3 |         8
#>    [4]         7        10         4 |         7
#>    [5]        11        15         5 |         6
#>    [6]        16        21         6 |         5
#>    [7]        22        28         7 |         4
#>    [8]        29        36         8 |         3
#>    [9]        37        45         9 |         2
#>   [10]        46        55        10 |         1

# must have names set
y <- IRanges::RleList(chr1 = x)
as_ranges(y)
#> GRanges object with 10 ranges and 1 metadata column:
#>        seqnames    ranges strand |     score
#>           <Rle> <IRanges>  <Rle> | <integer>
#>    [1]     chr1         1      * |        10
#>    [2]     chr1       2-3      * |         9
#>    [3]     chr1       4-6      * |         8
#>    [4]     chr1      7-10      * |         7
#>    [5]     chr1     11-15      * |         6
#>    [6]     chr1     16-21      * |         5
#>    [7]     chr1     22-28      * |         4
#>    [8]     chr1     29-36      * |         3
#>    [9]     chr1     37-45      * |         2
#>   [10]     chr1     46-55      * |         1
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
