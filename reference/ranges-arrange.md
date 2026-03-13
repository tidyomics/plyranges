# Sort a Ranges object

Sort a Ranges object

## Usage

``` r
# S3 method for class 'Ranges'
arrange(.data, ...)
```

## Arguments

- .data:

  A Ranges object.

- ...:

  Comma seperated list of variable names.

## Value

A sorted Ranges object

## Examples

``` r
rng <- as_iranges(data.frame(start = 1:10, width = 10:1))
rng <- mutate(rng, score = runif(10))
arrange(rng, score)
#> IRanges object with 10 ranges and 1 metadata column:
#>            start       end     width |     score
#>        <integer> <integer> <integer> | <numeric>
#>    [1]        10        10         1 | 0.0301457
#>    [2]         9        10         2 | 0.1252391
#>    [3]         3        10         8 | 0.2098168
#>    [4]         4        10         7 | 0.3581380
#>    [5]         7        10         4 | 0.3894393
#>    [6]         5        10         6 | 0.4482991
#>    [7]         8        10         3 | 0.5174597
#>    [8]         1        10        10 | 0.7334667
#>    [9]         6        10         5 | 0.9064264
#>   [10]         2        10         9 | 0.9069544
# you can also use dplyr::desc to arrange by descending order
```
