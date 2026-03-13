# Compute the number of distinct unique values in a vector or List

This is a wrapper to `length(unique(x))` or `lengths(unique(x))` if `x`
is a List object

## Usage

``` r
n_distinct(var)
```

## Arguments

- var:

  a vector of values

## Value

an integer vector

## Examples

``` r
x <- CharacterList(c("a", "b", "c", "a"),  "d")
n_distinct(x)
#> [1] 3 1
n_distinct(unlist(x))
#> [1] 4
```
