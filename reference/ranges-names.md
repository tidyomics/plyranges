# Tools for working with named Ranges

Tools for working with named Ranges

## Usage

``` r
remove_names(.data)

names_to_column(.data, var = "name")

id_to_column(.data, var = "id")
```

## Arguments

- .data:

  a Ranges object

- var:

  Name of column to use for names

## Value

Returns a Ranges object with empty names

## Details

The function `names_to_column()` and `id_to_column()` always places
`var` as the first column in `mcols(.data)`, shifting all other columns
to the left. The `id_to_column()` creates a column with sequential row
identifiers starting at 1, it will also remove any existing names.

## Examples

``` r
ir <- IRanges::IRanges(start = 1:3, width = 4, names = c("a", "b", "c"))
remove_names(ir)
#> IRanges object with 3 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]         1         4         4
#>   [2]         2         5         4
#>   [3]         3         6         4
ir_noname <- names_to_column(ir)
ir_noname
#> IRanges object with 3 ranges and 1 metadata column:
#>           start       end     width |        name
#>       <integer> <integer> <integer> | <character>
#>   [1]         1         4         4 |           a
#>   [2]         2         5         4 |           b
#>   [3]         3         6         4 |           c
ir_with_id <- id_to_column(ir)
ir_with_id
#> IRanges object with 3 ranges and 1 metadata column:
#>           start       end     width |        id
#>       <integer> <integer> <integer> | <integer>
#>   [1]         1         4         4 |         1
#>   [2]         2         5         4 |         2
#>   [3]         3         6         4 |         3
```
