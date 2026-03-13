# Vector-wise Range set-operations

Vector-wise Range set-operations

## Usage

``` r
intersect_ranges(x, y)

intersect_ranges_directed(x, y)

union_ranges(x, y)

union_ranges_directed(x, y)

setdiff_ranges(x, y)

setdiff_ranges_directed(x, y)

complement_ranges(x)

complement_ranges_directed(x)
```

## Arguments

- x, y:

  Two Ranges objects to compare.

## Value

A Ranges object

## Details

These are usual set-operations that act on the sets of the ranges
represented in x and y. By default these operations will ignore any
strand information. The directed versions of these functions will take
into account strand for GRanges objects.

## Examples

``` r
gr1 <- data.frame(seqnames = "chr1",
                  start = c(2,9),
                  end = c(7,9),
                  strand = c("+", "-")) %>%
               as_granges()
gr2 <- data.frame(seqnames = "chr1", start = 5, width = 5, strand = "-") %>%
         as_granges()

union_ranges(gr1, gr2)
#> GRanges object with 1 range and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1       2-9      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
union_ranges_directed(gr1, gr2)
#> GRanges object with 2 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1       2-7      +
#>   [2]     chr1       5-9      -
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

intersect_ranges(gr1, gr2)
#> GRanges object with 2 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1       5-7      *
#>   [2]     chr1         9      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
intersect_ranges_directed(gr1, gr2)
#> GRanges object with 1 range and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1         9      -
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

setdiff_ranges(gr1, gr2)
#> GRanges object with 1 range and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1       2-4      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
setdiff_ranges_directed(gr1, gr2)
#> GRanges object with 1 range and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1       2-7      +
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
# taking the complement of a ranges requires annotation information
gr1 <- set_genome_info(gr1, seqlengths = 100)
complement_ranges(gr1)
#> GRanges object with 3 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1         1      *
#>   [2]     chr1         8      *
#>   [3]     chr1    10-100      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
