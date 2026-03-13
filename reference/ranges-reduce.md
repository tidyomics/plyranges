# Reduce then aggregate a Ranges object

Reduce then aggregate a Ranges object

## Usage

``` r
reduce_ranges(.data, min.gapwidth = 1L, ...)

reduce_ranges_directed(.data, min.gapwidth = 1L, ...)
```

## Arguments

- .data:

  a Ranges object to reduce

- min.gapwidth:

  Ranges separated by a gap of at least min.gapwidth positions are not
  merged.

- ...:

  Name-value pairs of summary functions.

## Value

a Ranges object with the

## Examples

``` r
set.seed(10)
df <- data.frame(start = sample(1:10), 
                 width = 5,  
                 seqnames = "seq1",
                 strand = sample(c("+", "-", "*"), 10, replace = TRUE), 
                 gc = runif(10))
                 
rng <- as_granges(df)
rng %>% reduce_ranges()
#> GRanges object with 1 range and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     seq1      1-14      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
rng %>% reduce_ranges(gc = mean(gc))
#> GRanges object with 1 range and 1 metadata column:
#>       seqnames    ranges strand |        gc
#>          <Rle> <IRanges>  <Rle> | <numeric>
#>   [1]     seq1      1-14      * |  0.583196
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
rng %>% reduce_ranges_directed(gc = mean(gc))
#> GRanges object with 4 ranges and 1 metadata column:
#>       seqnames    ranges strand |        gc
#>          <Rle> <IRanges>  <Rle> | <numeric>
#>   [1]     seq1       1-5      + |  0.500503
#>   [2]     seq1      9-13      + |  0.535597
#>   [3]     seq1      3-12      - |  0.558135
#>   [4]     seq1      2-14      * |  0.640830
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
rng %>% reduce_ranges_directed(gc = mean(gc), min.gapwidth = 10)
#> GRanges object with 3 ranges and 1 metadata column:
#>       seqnames    ranges strand |        gc
#>          <Rle> <IRanges>  <Rle> | <numeric>
#>   [1]     seq1      1-13      + |  0.518050
#>   [2]     seq1      3-12      - |  0.558135
#>   [3]     seq1      2-14      * |  0.640830
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

x <- data.frame(start = c(11:13, 2, 7:6), 
               width=3, 
               id=sample(letters[1:3], 6, replace = TRUE),
               score= sample(1:6))
x <- as_iranges(x)
x %>% reduce_ranges()
#> IRanges object with 3 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]         2         4         3
#>   [2]         6         9         4
#>   [3]        11        15         5
x %>% reduce_ranges(score = sum(score))
#> IRanges object with 3 ranges and 1 metadata column:
#>           start       end     width |     score
#>       <integer> <integer> <integer> | <integer>
#>   [1]         2         4         3 |         5
#>   [2]         6         9         4 |         5
#>   [3]        11        15         5 |        11
x %>% group_by(id) %>% reduce_ranges(score = sum(score))
#> IRanges object with 4 ranges and 2 metadata columns:
#>           start       end     width |          id     score
#>       <integer> <integer> <integer> | <character> <integer>
#>   [1]        11        13         3 |           c         2
#>   [2]         2         4         3 |           b         5
#>   [3]        12        15         4 |           b         9
#>   [4]         6         9         4 |           a         5
```
