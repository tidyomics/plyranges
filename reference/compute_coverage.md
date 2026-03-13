# Compute coverage over a Ranges object

Compute coverage over a Ranges object

## Usage

``` r
compute_coverage(x, shift, width, weight, ...)
```

## Arguments

- x:

  a `Ranges` object

- shift:

  shift how much should each range in x be shifted by? (default = 0L)

- width:

  width how long should the returned coverage score be? This must be
  either a positive integer or NULL (default = NULL)

- weight:

  weight how much weight should be assigned to each range? Either an
  integer or numeric vector or a column in x. (default = 1L)

- ...:

  other optional parameters to pass to coverage

## Value

An expanded Ranges object with a score column corresponding to the
coverage value over that interval. Note that compute_coverage drops
metadata associated with the orginal ranges.

## See also

`IRanges::`[`coverage()`](https://rdrr.io/pkg/IRanges/man/coverage-methods.html),
`GenomicRanges::`[`coverage()`](https://rdrr.io/pkg/GenomicRanges/man/coverage-methods.html)

## Examples

``` r
rng <- as_iranges(data.frame(start = 1:10, width = 5))
compute_coverage(rng)
#> IRanges object with 9 ranges and 1 metadata column:
#>           start       end     width |     score
#>       <integer> <integer> <integer> | <integer>
#>   [1]         1         1         1 |         1
#>   [2]         2         2         1 |         2
#>   [3]         3         3         1 |         3
#>   [4]         4         4         1 |         4
#>   [5]         5        10         6 |         5
#>   [6]        11        11         1 |         4
#>   [7]        12        12         1 |         3
#>   [8]        13        13         1 |         2
#>   [9]        14        14         1 |         1
compute_coverage(rng, shift = 14L)
#> IRanges object with 10 ranges and 1 metadata column:
#>            start       end     width |     score
#>        <integer> <integer> <integer> | <integer>
#>    [1]         1        14        14 |         0
#>    [2]        15        15         1 |         1
#>    [3]        16        16         1 |         2
#>    [4]        17        17         1 |         3
#>    [5]        18        18         1 |         4
#>    [6]        19        24         6 |         5
#>    [7]        25        25         1 |         4
#>    [8]        26        26         1 |         3
#>    [9]        27        27         1 |         2
#>   [10]        28        28         1 |         1
compute_coverage(rng, width = 10L)
#> IRanges object with 5 ranges and 1 metadata column:
#>           start       end     width |     score
#>       <integer> <integer> <integer> | <integer>
#>   [1]         1         1         1 |         1
#>   [2]         2         2         1 |         2
#>   [3]         3         3         1 |         3
#>   [4]         4         4         1 |         4
#>   [5]         5        10         6 |         5
```
