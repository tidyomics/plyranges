# Add distance to nearest neighbours between two Ranges objects

Appends distance to nearest subject range to query ranges similar to
setting `distance` in `join_nearest_`. Distance is set to `NA` for
features with no nearest feature by the selected nearest metric.

## Usage

``` r
add_nearest_distance(x, y = x, name = "distance")

add_nearest_distance_left(x, y = x, name = "distance")

add_nearest_distance_right(x, y = x, name = "distance")

add_nearest_distance_upstream(x, y = x, name = "distance")

add_nearest_distance_downstream(x, y = x, name = "distance")
```

## Arguments

- x:

  The query ranges

- y:

  the subject ranges within which the nearest ranges are found. If
  missing, query ranges are used as the subject.

- name:

  column name to create containing distance values

## Value

ranges in `x` with additional column containing the distance to the
nearest range in `y`.

## Details

By default `add_nearest_distance` will find arbitrary nearest neighbours
in either direction and ignore any strand information. The
`add_nearest_distance_left` and `add_nearest_distance_right` methods
will find arbitrary nearest neighbour ranges on x that are left/right of
those on y and ignore any strand information.

The `add_nearest_distance_upstream` method will find arbitrary nearest
neighbour ranges on x that are upstream of those on y. This takes into
account strandedness of the ranges. On the positive strand nearest
upstream will be on the left and on the negative strand nearest upstream
will be on the right.

The `add_nearest_distance_downstream` method will find arbitrary nearest
neighbour ranges on x that are upstream of those on y. This takes into
account strandedness of the ranges. On the positive strand nearest
downstream will be on the right and on the negative strand nearest
upstream will be on the left.

## See also

[`join_nearest`](https://tidyomics.github.io/plyranges/reference/ranges-nearest.md)

## Examples

``` r
query <- data.frame(start = c(5,10, 15,20),
                   width = 5,
                   gc = runif(4)) %>%
             as_iranges()
subject <- data.frame(start = c(2:6, 24),
                      width = 3:8,
                      label = letters[1:6]) %>%
             as_iranges()
             
add_nearest_distance(query, subject)
#> IRanges object with 4 ranges and 2 metadata columns:
#>           start       end     width |        gc  distance
#>       <integer> <integer> <integer> | <numeric> <integer>
#>   [1]         5         9         5 | 0.0807501         0
#>   [2]        10        14         5 | 0.8343330         0
#>   [3]        15        19         5 | 0.6007609         2
#>   [4]        20        24         5 | 0.1572084         0
add_nearest_distance_left(query, subject)
#> IRanges object with 4 ranges and 2 metadata columns:
#>           start       end     width |        gc  distance
#>       <integer> <integer> <integer> | <numeric> <integer>
#>   [1]         5         9         5 | 0.0807501         0
#>   [2]        10        14         5 | 0.8343330         0
#>   [3]        15        19         5 | 0.6007609         2
#>   [4]        20        24         5 | 0.1572084      <NA>
add_nearest_distance_left(query)
#> IRanges object with 4 ranges and 2 metadata columns:
#>           start       end     width |        gc  distance
#>       <integer> <integer> <integer> | <numeric> <integer>
#>   [1]         5         9         5 | 0.0807501      <NA>
#>   [2]        10        14         5 | 0.8343330         0
#>   [3]        15        19         5 | 0.6007609         0
#>   [4]        20        24         5 | 0.1572084         0
```
