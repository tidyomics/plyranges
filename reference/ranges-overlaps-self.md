# Find overlaps within a Ranges object

Find overlaps within a Ranges object

## Usage

``` r
join_overlap_self(x, maxgap, minoverlap)

join_overlap_self_within(x, maxgap, minoverlap)

join_overlap_self_directed(x, maxgap, minoverlap)

join_overlap_self_within_directed(x, maxgap, minoverlap)
```

## Arguments

- x:

  A Ranges object

- maxgap, minoverlap:

  The maximimum gap between intervals as an integer greater than or
  equal to zero. The minimum amount of overlap between intervals as an
  integer greater than zero, accounting for the maximum gap.

## Value

a Ranges object

## Details

Self overlaps find any overlaps (or overlaps within or overlaps
directed) between a ranges object and itself.

## See also

[`find_overlaps()`](https://tidyomics.github.io/plyranges/reference/ranges-overlaps.md),
[`join_overlap_inner()`](https://tidyomics.github.io/plyranges/reference/overlap-joins.md)

## Examples

``` r
query <- data.frame(start = c(5,10, 15,20), width = 5, gc = runif(4)) %>%
             as_iranges()

join_overlap_self(query)
#> IRanges object with 4 ranges and 2 metadata columns:
#>           start       end     width |        gc gc.overlap
#>       <integer> <integer> <integer> | <numeric>  <numeric>
#>   [1]         5         9         5 | 0.0161834  0.0161834
#>   [2]        10        14         5 | 0.4721299  0.4721299
#>   [3]        15        19         5 | 0.0423897  0.0423897
#>   [4]        20        24         5 | 0.4631544  0.4631544

# -- GRanges objects, strand is ignored by default
query  <- data.frame(seqnames = "chr1",
               start = c(11,101),
               end = c(21, 200),
               name = c("a1", "a2"),
               strand = c("+", "-"),
               score = c(1,2)) %>%
           as_granges()

# ignores strandedness
join_overlap_self(query)
#> GRanges object with 2 ranges and 4 metadata columns:
#>       seqnames    ranges strand |        name     score name.overlap
#>          <Rle> <IRanges>  <Rle> | <character> <numeric>  <character>
#>   [1]     chr1     11-21      + |          a1         1           a1
#>   [2]     chr1   101-200      - |          a2         2           a2
#>       score.overlap
#>           <numeric>
#>   [1]             1
#>   [2]             2
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
join_overlap_self_within(query)
#> GRanges object with 2 ranges and 4 metadata columns:
#>       seqnames    ranges strand |        name     score name.overlap
#>          <Rle> <IRanges>  <Rle> | <character> <numeric>  <character>
#>   [1]     chr1     11-21      + |          a1         1           a1
#>   [2]     chr1   101-200      - |          a2         2           a2
#>       score.overlap
#>           <numeric>
#>   [1]             1
#>   [2]             2
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
# adding directed prefix includes strand
join_overlap_self_directed(query)
#> GRanges object with 2 ranges and 4 metadata columns:
#>       seqnames    ranges strand |        name     score name.overlap
#>          <Rle> <IRanges>  <Rle> | <character> <numeric>  <character>
#>   [1]     chr1     11-21      + |          a1         1           a1
#>   [2]     chr1   101-200      - |          a2         2           a2
#>       score.overlap
#>           <numeric>
#>   [1]             1
#>   [2]             2
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

```
