# Compute the number of ranges in each group.

This function should only be used within `summarise()`, `mutate()` and
[`filter()`](https://rdrr.io/r/stats/filter.html).

## Usage

``` r
n()
```

## Value

`n()` will only be evaluated inside a function call, where it returns an
integer.

## Examples

``` r
ir <- as_iranges(
                 data.frame(start = 1:10,
                            width = 5,
                            name = c(rep("a", 5), rep("b", 3), rep("c", 2))
                            )
                )
by_names <- group_by(ir, name)
summarise(by_names, n = n())
#> DataFrame with 3 rows and 2 columns
#>          name         n
#>   <character> <integer>
#> 1           a         5
#> 2           b         3
#> 3           c         2
mutate(by_names, n = n())
#> IRanges object with 10 ranges and 2 metadata columns:
#> Groups: name [3]
#>            start       end     width |        name         n
#>        <integer> <integer> <integer> | <character> <integer>
#>    [1]         1         5         5 |           a         5
#>    [2]         2         6         5 |           a         5
#>    [3]         3         7         5 |           a         5
#>    [4]         4         8         5 |           a         5
#>    [5]         5         9         5 |           a         5
#>    [6]         6        10         5 |           b         3
#>    [7]         7        11         5 |           b         3
#>    [8]         8        12         5 |           b         3
#>    [9]         9        13         5 |           c         2
#>   [10]        10        14         5 |           c         2
filter(by_names, n() >= 3)
#> IRanges object with 8 ranges and 1 metadata column:
#> Groups: name [2]
#>           start       end     width |        name
#>       <integer> <integer> <integer> | <character>
#>   [1]         1         5         5 |           a
#>   [2]         2         6         5 |           a
#>   [3]         3         7         5 |           a
#>   [4]         4         8         5 |           a
#>   [5]         5         9         5 |           a
#>   [6]         6        10         5 |           b
#>   [7]         7        11         5 |           b
#>   [8]         8        12         5 |           b
```
