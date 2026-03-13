# Disjoin then aggregate a Ranges object

Disjoin then aggregate a Ranges object

## Usage

``` r
disjoin_ranges(.data, ...)

disjoin_ranges_directed(.data, ...)
```

## Arguments

- .data:

  a Ranges object to disjoin

- ...:

  Name-value pairs of summary functions.

## Value

a Ranges object that is now disjoint (no bases overlap).

## Examples

``` r
df <- data.frame(start = 1:10, width = 5,  seqnames = "seq1",
strand = sample(c("+", "-", "*"), 10, replace = TRUE), gc = runif(10))
rng <- as_granges(df)
rng %>% disjoin_ranges()
#> GRanges object with 14 ranges and 0 metadata columns:
#>        seqnames    ranges strand
#>           <Rle> <IRanges>  <Rle>
#>    [1]     seq1         1      *
#>    [2]     seq1         2      *
#>    [3]     seq1         3      *
#>    [4]     seq1         4      *
#>    [5]     seq1         5      *
#>    ...      ...       ...    ...
#>   [10]     seq1        10      *
#>   [11]     seq1        11      *
#>   [12]     seq1        12      *
#>   [13]     seq1        13      *
#>   [14]     seq1        14      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
rng %>% disjoin_ranges(gc = mean(gc))
#> GRanges object with 14 ranges and 1 metadata column:
#>        seqnames    ranges strand |        gc
#>           <Rle> <IRanges>  <Rle> | <numeric>
#>    [1]     seq1         1      * |  0.726227
#>    [2]     seq1         2      * |  0.463501
#>    [3]     seq1         3      * |  0.589068
#>    [4]     seq1         4      * |  0.540966
#>    [5]     seq1         5      * |  0.511314
#>    ...      ...       ...    ... .       ...
#>   [10]     seq1        10      * |  0.486433
#>   [11]     seq1        11      * |  0.489902
#>   [12]     seq1        12      * |  0.458650
#>   [13]     seq1        13      * |  0.511781
#>   [14]     seq1        14      * |  0.995088
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
rng %>% disjoin_ranges_directed(gc = mean(gc))
#> GRanges object with 14 ranges and 1 metadata column:
#>        seqnames    ranges strand |        gc
#>           <Rle> <IRanges>  <Rle> | <numeric>
#>    [1]     seq1         1      + |  0.726227
#>    [2]     seq1         2      + |  0.463501
#>    [3]     seq1         3      + |  0.589068
#>    [4]     seq1       4-5      + |  0.540966
#>    [5]     seq1         6      + |  0.479213
#>    ...      ...       ...    ... .       ...
#>   [10]     seq1       5-7      - |  0.392702
#>   [11]     seq1       8-9      - |  0.372545
#>   [12]     seq1     10-12      - |  0.673738
#>   [13]     seq1     13-14      - |  0.995088
#>   [14]     seq1      6-10      * |  0.472557
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
