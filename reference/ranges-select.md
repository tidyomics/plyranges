# Select metadata columns of the Ranges object by name or position

Select metadata columns of the Ranges object by name or position

## Usage

``` r
# S3 method for class 'Ranges'
select(.data, ..., .drop_ranges = FALSE)
```

## Arguments

- .data:

  a `Ranges` object

- ...:

  One or more metadata column names.

- .drop_ranges:

  If TRUE select will always return a tibble. In this case, you may
  select columns that form the core part of the Ranges object.

## Value

a Ranges object or a tibble

## Details

Note that by default select only acts on the metadata columns (and will
therefore return a Ranges object) if a core component of a Ranges is
dropped or selected without the other required components (this includes
the seqnames, strand, start, end, width names), then select will throw
an error unless .drop_ranges is set to TRUE.

## See also

[`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html)

## Examples

``` r
df <- data.frame(start = 1:10, width = 5,  seqnames = "seq1",
strand = sample(c("+", "-", "*"), 10, replace = TRUE), gc = runif(10), counts = rpois(10, 2))
rng <- as_granges(df)
select(rng, -gc)
#> GRanges object with 10 ranges and 1 metadata column:
#>        seqnames    ranges strand |    counts
#>           <Rle> <IRanges>  <Rle> | <integer>
#>    [1]     seq1       1-5      * |         2
#>    [2]     seq1       2-6      * |         1
#>    [3]     seq1       3-7      + |         5
#>    [4]     seq1       4-8      - |         1
#>    [5]     seq1       5-9      * |         1
#>    [6]     seq1      6-10      - |         4
#>    [7]     seq1      7-11      - |         2
#>    [8]     seq1      8-12      + |         1
#>    [9]     seq1      9-13      + |         2
#>   [10]     seq1     10-14      + |         1
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
select(rng, gc)
#> GRanges object with 10 ranges and 1 metadata column:
#>        seqnames    ranges strand |        gc
#>           <Rle> <IRanges>  <Rle> | <numeric>
#>    [1]     seq1       1-5      * | 0.7611217
#>    [2]     seq1       2-6      * | 0.5733564
#>    [3]     seq1       3-7      + | 0.4475080
#>    [4]     seq1       4-8      - | 0.0838020
#>    [5]     seq1       5-9      * | 0.2191385
#>    [6]     seq1      6-10      - | 0.0755703
#>    [7]     seq1      7-11      - | 0.5344268
#>    [8]     seq1      8-12      + | 0.6413566
#>    [9]     seq1      9-13      + | 0.5257393
#>   [10]     seq1     10-14      + | 0.0392814
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
select(rng, counts, gc)
#> GRanges object with 10 ranges and 2 metadata columns:
#>        seqnames    ranges strand |    counts        gc
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric>
#>    [1]     seq1       1-5      * |         2 0.7611217
#>    [2]     seq1       2-6      * |         1 0.5733564
#>    [3]     seq1       3-7      + |         5 0.4475080
#>    [4]     seq1       4-8      - |         1 0.0838020
#>    [5]     seq1       5-9      * |         1 0.2191385
#>    [6]     seq1      6-10      - |         4 0.0755703
#>    [7]     seq1      7-11      - |         2 0.5344268
#>    [8]     seq1      8-12      + |         1 0.6413566
#>    [9]     seq1      9-13      + |         2 0.5257393
#>   [10]     seq1     10-14      + |         1 0.0392814
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
select(rng, 2:1)
#> GRanges object with 10 ranges and 2 metadata columns:
#>        seqnames    ranges strand |    counts        gc
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric>
#>    [1]     seq1       1-5      * |         2 0.7611217
#>    [2]     seq1       2-6      * |         1 0.5733564
#>    [3]     seq1       3-7      + |         5 0.4475080
#>    [4]     seq1       4-8      - |         1 0.0838020
#>    [5]     seq1       5-9      * |         1 0.2191385
#>    [6]     seq1      6-10      - |         4 0.0755703
#>    [7]     seq1      7-11      - |         2 0.5344268
#>    [8]     seq1      8-12      + |         1 0.6413566
#>    [9]     seq1      9-13      + |         2 0.5257393
#>   [10]     seq1     10-14      + |         1 0.0392814
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
select(rng, seqnames, strand, .drop_ranges = TRUE)
#> DataFrame with 10 rows and 2 columns
#>    seqnames strand
#>       <Rle>  <Rle>
#> 1      seq1      *
#> 2      seq1      *
#> 3      seq1      +
#> 4      seq1      -
#> 5      seq1      *
#> 6      seq1      -
#> 7      seq1      -
#> 8      seq1      +
#> 9      seq1      +
#> 10     seq1      +
```
