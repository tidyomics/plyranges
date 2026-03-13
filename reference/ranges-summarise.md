# Reduce multiple values in a Ranges down to a single value

Reduce multiple values in a Ranges down to a single value

## Usage

``` r
# S3 method for class 'Ranges'
summarise(.data, ...)
```

## Arguments

- .data:

  a Ranges object

- ...:

  Name-value pairs of summary functions. The name will be the name of
  the variable in the result. The value should be an expression that
  will return a value that has length one or length equal to the number
  of groups.

## Value

A
`S4Vectors::`[`DataFrame()`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)

## Details

Creates one or more variables as a
`S4Vectors::`[`DataFrame()`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
from the input Ranges object. If the ranges object is grouped, there
will be a row for each group. Because grouping may remove whether a
Ranges object is valid, a DataFrame is always returned.

## Examples

``` r
df <- data.frame(start = 1:10, width = 5,  seqnames = "seq1",
strand = sample(c("+", "-", "*"), 10, replace = TRUE), gc = runif(10))
rng <- as_granges(df)
rng %>% summarise(gc = mean(gc))
#> DataFrame with 1 row and 1 column
#>          gc
#>   <numeric>
#> 1  0.454874
rng %>% group_by(strand) %>% summarise(gc = mean(gc))
#> DataFrame with 3 rows and 2 columns
#>   strand        gc
#>    <Rle> <numeric>
#> 1      +  0.648149
#> 2      -  0.475607
#> 3      *  0.337504
```
