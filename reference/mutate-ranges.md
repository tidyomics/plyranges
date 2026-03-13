# Modify a Ranges object

Modify a Ranges object

## Usage

``` r
# S3 method for class 'Ranges'
mutate(.data, ...)
```

## Arguments

- .data:

  a `Ranges` object

- ...:

  Pairs of name-value expressions. The name-value pairs can either
  create new metadata columns or modify existing ones.

## Value

a Ranges object

## Examples

``` r
df <- data.frame(start = 1:10,
                 width = 5,
                 seqnames = "seq1",
                 strand = sample(c("+", "-", "*"), 10, replace = TRUE),
                 gc = runif(10))
rng <- as_granges(df)

# mutate adds new columns
rng %>%
    mutate(avg_gc = mean(gc), row_id = 1:n())
#> GRanges object with 10 ranges and 3 metadata columns:
#>        seqnames    ranges strand |        gc    avg_gc    row_id
#>           <Rle> <IRanges>  <Rle> | <numeric> <numeric> <integer>
#>    [1]     seq1       1-5      - |  0.884227  0.340702         1
#>    [2]     seq1       2-6      * |  0.207714  0.340702         2
#>    [3]     seq1       3-7      - |  0.307086  0.340702         3
#>    [4]     seq1       4-8      + |  0.330530  0.340702         4
#>    [5]     seq1       5-9      * |  0.198679  0.340702         5
#>    [6]     seq1      6-10      * |  0.235694  0.340702         6
#>    [7]     seq1      7-11      + |  0.274887  0.340702         7
#>    [8]     seq1      8-12      * |  0.591321  0.340702         8
#>    [9]     seq1      9-13      - |  0.253391  0.340702         9
#>   [10]     seq1     10-14      - |  0.123487  0.340702        10
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
# can also compute on newly created columns
rng %>%
    mutate(score = gc * width, score2 = score + 1)
#> GRanges object with 10 ranges and 3 metadata columns:
#>        seqnames    ranges strand |        gc     score    score2
#>           <Rle> <IRanges>  <Rle> | <numeric> <numeric> <numeric>
#>    [1]     seq1       1-5      - |  0.884227  4.421135   5.42114
#>    [2]     seq1       2-6      * |  0.207714  1.038569   2.03857
#>    [3]     seq1       3-7      - |  0.307086  1.535429   2.53543
#>    [4]     seq1       4-8      + |  0.330530  1.652649   2.65265
#>    [5]     seq1       5-9      * |  0.198679  0.993395   1.99340
#>    [6]     seq1      6-10      * |  0.235694  1.178472   2.17847
#>    [7]     seq1      7-11      + |  0.274887  1.374433   2.37443
#>    [8]     seq1      8-12      * |  0.591321  2.956605   3.95661
#>    [9]     seq1      9-13      - |  0.253391  1.266953   2.26695
#>   [10]     seq1     10-14      - |  0.123487  0.617436   1.61744
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
# group by partitions the data and computes within each group
rng %>%
    group_by(strand) %>%
    mutate(avg_gc = mean(gc), row_id = 1:n())
#> GRanges object with 10 ranges and 3 metadata columns:
#> Groups: strand [3]
#>        seqnames    ranges strand |        gc    avg_gc    row_id
#>           <Rle> <IRanges>  <Rle> | <numeric> <numeric> <integer>
#>    [1]     seq1       1-5      - |  0.884227  0.392048         1
#>    [2]     seq1       2-6      * |  0.207714  0.308352         1
#>    [3]     seq1       3-7      - |  0.307086  0.392048         2
#>    [4]     seq1       4-8      + |  0.330530  0.302708         1
#>    [5]     seq1       5-9      * |  0.198679  0.308352         2
#>    [6]     seq1      6-10      * |  0.235694  0.308352         3
#>    [7]     seq1      7-11      + |  0.274887  0.302708         2
#>    [8]     seq1      8-12      * |  0.591321  0.308352         4
#>    [9]     seq1      9-13      - |  0.253391  0.392048         3
#>   [10]     seq1     10-14      - |  0.123487  0.392048         4
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

# mutate can be used in conjuction with anchoring to resize ranges
rng %>%
    mutate(width = 10)
#> GRanges object with 10 ranges and 1 metadata column:
#>        seqnames    ranges strand |        gc
#>           <Rle> <IRanges>  <Rle> | <numeric>
#>    [1]     seq1      1-10      - |  0.884227
#>    [2]     seq1      2-11      * |  0.207714
#>    [3]     seq1      3-12      - |  0.307086
#>    [4]     seq1      4-13      + |  0.330530
#>    [5]     seq1      5-14      * |  0.198679
#>    [6]     seq1      6-15      * |  0.235694
#>    [7]     seq1      7-16      + |  0.274887
#>    [8]     seq1      8-17      * |  0.591321
#>    [9]     seq1      9-18      - |  0.253391
#>   [10]     seq1     10-19      - |  0.123487
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
# by default width modfication fixes by start
rng %>%
    anchor_start() %>%
    mutate(width = 10)
#> GRanges object with 10 ranges and 1 metadata column:
#>        seqnames    ranges strand |        gc
#>           <Rle> <IRanges>  <Rle> | <numeric>
#>    [1]     seq1      1-10      - |  0.884227
#>    [2]     seq1      2-11      * |  0.207714
#>    [3]     seq1      3-12      - |  0.307086
#>    [4]     seq1      4-13      + |  0.330530
#>    [5]     seq1      5-14      * |  0.198679
#>    [6]     seq1      6-15      * |  0.235694
#>    [7]     seq1      7-16      + |  0.274887
#>    [8]     seq1      8-17      * |  0.591321
#>    [9]     seq1      9-18      - |  0.253391
#>   [10]     seq1     10-19      - |  0.123487
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
# fix by end or midpoint
rng %>%
    anchor_end() %>%
    mutate(width = width + 1)
#> GRanges object with 10 ranges and 1 metadata column:
#>        seqnames    ranges strand |        gc
#>           <Rle> <IRanges>  <Rle> | <numeric>
#>    [1]     seq1       0-5      - |  0.884227
#>    [2]     seq1       1-6      * |  0.207714
#>    [3]     seq1       2-7      - |  0.307086
#>    [4]     seq1       3-8      + |  0.330530
#>    [5]     seq1       4-9      * |  0.198679
#>    [6]     seq1      5-10      * |  0.235694
#>    [7]     seq1      6-11      + |  0.274887
#>    [8]     seq1      7-12      * |  0.591321
#>    [9]     seq1      8-13      - |  0.253391
#>   [10]     seq1      9-14      - |  0.123487
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
rng %>%
    anchor_center() %>%
    mutate(width = width + 1)
#> GRanges object with 10 ranges and 1 metadata column:
#>        seqnames    ranges strand |        gc
#>           <Rle> <IRanges>  <Rle> | <numeric>
#>    [1]     seq1       0-5      - |  0.884227
#>    [2]     seq1       1-6      * |  0.207714
#>    [3]     seq1       2-7      - |  0.307086
#>    [4]     seq1       3-8      + |  0.330530
#>    [5]     seq1       4-9      * |  0.198679
#>    [6]     seq1      5-10      * |  0.235694
#>    [7]     seq1      6-11      + |  0.274887
#>    [8]     seq1      7-12      * |  0.591321
#>    [9]     seq1      8-13      - |  0.253391
#>   [10]     seq1      9-14      - |  0.123487
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
# anchoring by strand
rng %>%
    anchor_3p() %>%
    mutate(width = width * 2)
#> GRanges object with 10 ranges and 1 metadata column:
#>        seqnames    ranges strand |        gc
#>           <Rle> <IRanges>  <Rle> | <numeric>
#>    [1]     seq1      1-10      - |  0.884227
#>    [2]     seq1      -3-6      * |  0.207714
#>    [3]     seq1      3-12      - |  0.307086
#>    [4]     seq1      -1-8      + |  0.330530
#>    [5]     seq1       0-9      * |  0.198679
#>    [6]     seq1      1-10      * |  0.235694
#>    [7]     seq1      2-11      + |  0.274887
#>    [8]     seq1      3-12      * |  0.591321
#>    [9]     seq1      9-18      - |  0.253391
#>   [10]     seq1     10-19      - |  0.123487
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
rng %>%
    anchor_5p() %>%
    mutate(width = width * 2)
#> GRanges object with 10 ranges and 1 metadata column:
#>        seqnames    ranges strand |        gc
#>           <Rle> <IRanges>  <Rle> | <numeric>
#>    [1]     seq1      -4-5      - |  0.884227
#>    [2]     seq1      2-11      * |  0.207714
#>    [3]     seq1      -2-7      - |  0.307086
#>    [4]     seq1      4-13      + |  0.330530
#>    [5]     seq1      5-14      * |  0.198679
#>    [6]     seq1      6-15      * |  0.235694
#>    [7]     seq1      7-16      + |  0.274887
#>    [8]     seq1      8-17      * |  0.591321
#>    [9]     seq1      4-13      - |  0.253391
#>   [10]     seq1      5-14      - |  0.123487
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
