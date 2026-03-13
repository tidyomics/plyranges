# Slide or tile over a Ranges object

Slide or tile over a Ranges object

## Usage

``` r
tile_ranges(x, width)

slide_ranges(x, width, step)
```

## Arguments

- x:

  a Ranges object

- width:

  the maximum width of each window/tile (integer vector of length 1)

- step:

  the distance between start position of each sliding window (integer
  vector of length 1)

## Value

a Ranges object

## Details

The `tile_ranges()` function paritions a Ranges object `x` by the given
the `width` over all ranges in `x`, truncated by the sequence end. The
`slide_ranges()` function makes sliding windows within each range of `x`
of size `width` and sliding by `step`. Both `slide_ranges()` and
`tile_ranges()` return a new Ranges object with a metadata column called
"partition" which contains the index of the input range `x` that a
parition belongs to.

## See also

`GenomicRanges::`[`tile()`](https://rdrr.io/pkg/GenomicRanges/man/tile-methods.html)

## Examples

``` r
 
gr <- data.frame(seqnames = c("chr1", rep("chr2", 3), rep("chr1", 2), rep("chr3", 4)),
                 start = 1:10,
                 end = 11,
                 strand = c("-", rep("+", 2), rep("*", 2), rep("+", 3), rep("-", 2))) %>%
      as_granges() %>%
      set_genome_info(seqlengths = c(11,12,13))

# partition ranges into subranges of width 2, odd width ranges
# will have one subrange of width 1              
tile_ranges(gr, width = 2)
#> GRanges object with 35 ranges and 1 metadata column:
#>        seqnames    ranges strand | partition
#>           <Rle> <IRanges>  <Rle> | <integer>
#>    [1]     chr1         1      - |         1
#>    [2]     chr1       2-3      - |         1
#>    [3]     chr1       4-5      - |         1
#>    [4]     chr1       6-7      - |         1
#>    [5]     chr1       8-9      - |         1
#>    ...      ...       ...    ... .       ...
#>   [31]     chr3       8-9      + |         8
#>   [32]     chr3     10-11      + |         8
#>   [33]     chr3         9      - |         9
#>   [34]     chr3     10-11      - |         9
#>   [35]     chr3     10-11      - |        10
#>   -------
#>   seqinfo: 3 sequences from an unspecified genome; no seqlengths

# make sliding windows of width 3, moving window with step size of 2
slide_ranges(gr, width = 3, step = 2)
#> GRanges object with 30 ranges and 1 metadata column:
#>        seqnames    ranges strand | partition
#>           <Rle> <IRanges>  <Rle> | <integer>
#>    [1]     chr1       1-3      - |         1
#>    [2]     chr1       3-5      - |         1
#>    [3]     chr1       5-7      - |         1
#>    [4]     chr1       7-9      - |         1
#>    [5]     chr1      9-11      - |         1
#>    ...      ...       ...    ... .       ...
#>   [26]     chr3      9-11      + |         7
#>   [27]     chr3      8-10      + |         8
#>   [28]     chr3     10-11      + |         8
#>   [29]     chr3      9-11      - |         9
#>   [30]     chr3     10-11      - |        10
#>   -------
#>   seqinfo: 3 sequences from an unspecified genome
```
