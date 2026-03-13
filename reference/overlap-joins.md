# Join by overlapping Ranges

Join by overlapping Ranges

## Usage

``` r
join_overlap_intersect(x, y, maxgap, minoverlap, suffix = c(".x", ".y"))

join_overlap_intersect_within(x, y, maxgap, minoverlap, suffix = c(".x", ".y"))

join_overlap_intersect_directed(
  x,
  y,
  maxgap,
  minoverlap,
  suffix = c(".x", ".y")
)

join_overlap_intersect_within_directed(
  x,
  y,
  maxgap,
  minoverlap,
  suffix = c(".x", ".y")
)

join_overlap_inner(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y"))

join_overlap_inner_within(
  x,
  y,
  maxgap = -1L,
  minoverlap = 0L,
  suffix = c(".x", ".y")
)

join_overlap_inner_directed(
  x,
  y,
  maxgap = -1L,
  minoverlap = 0L,
  suffix = c(".x", ".y")
)

join_overlap_inner_within_directed(
  x,
  y,
  maxgap = -1L,
  minoverlap = 0L,
  suffix = c(".x", ".y")
)

join_overlap_left(x, y, maxgap, minoverlap, suffix = c(".x", ".y"))

join_overlap_left_within(x, y, maxgap, minoverlap, suffix = c(".x", ".y"))

join_overlap_left_directed(x, y, maxgap, minoverlap, suffix = c(".x", ".y"))

join_overlap_left_within_directed(
  x,
  y,
  maxgap,
  minoverlap,
  suffix = c(".x", ".y")
)
```

## Arguments

- x, y:

  Objects representing ranges

- maxgap, minoverlap:

  The maximimum gap between intervals as an integer greater than or
  equal to zero. The minimum amount of overlap between intervals as an
  integer greater than zero, accounting for the maximum gap.

- suffix:

  Character to vectors to append to common columns in x and y (default =
  `c(".x", ".y")`).

## Value

a GRanges object

## Details

The function `join_overlap_intersect()` finds the genomic intervals that
are the overlapping ranges between x and y and returns a new ranges
object with metadata columns from x and y.

The function `join_overlap_inner()` is equivalent to
[`find_overlaps()`](https://tidyomics.github.io/plyranges/reference/ranges-overlaps.md).

The function `join_overlap_left()` performs a left outer join between x
and y. It returns all ranges in x that overlap or do not overlap ranges
in y plus metadata columns common to both. If there is no overlapping
range the metadata column will contain a missing value.

The function
[`join_overlap_self()`](https://tidyomics.github.io/plyranges/reference/ranges-overlaps-self.md)
find all overlaps between a ranges object x and itself.

All of these functions have two suffixes that modify their behavior. The
`within` suffix, returns only ranges in x that are completely overlapped
within in y. The `directed` suffix accounts for the strandedness of the
ranges when performing overlaps.

## See also

[`join_overlap_self()`](https://tidyomics.github.io/plyranges/reference/ranges-overlaps-self.md),
`join_overlap_left()`,
[`find_overlaps()`](https://tidyomics.github.io/plyranges/reference/ranges-overlaps.md)

## Examples

``` r
x <- as_iranges(data.frame(start = c(11, 101), end = c(21, 201)))
y <- as_iranges(data.frame(start = c(10, 20, 50, 100, 1),
                           end = c(19, 21, 105, 202, 5)))

# self
join_overlap_self(y)
#> IRanges object with 7 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]        10        19        10
#>   [2]        20        21         2
#>   [3]        50       105        56
#>   [4]        50       105        56
#>   [5]       100       202       103
#>   [6]       100       202       103
#>   [7]         1         5         5

# intersect takes common interval
join_overlap_intersect(x,y)
#> IRanges object with 4 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]        11        19         9
#>   [2]        20        21         2
#>   [3]       101       105         5
#>   [4]       101       201       101

# within
join_overlap_intersect_within(x,y)
#> IRanges object with 1 range and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]       101       201       101

# left, and inner join, it's often useful having an id column here
y <- y %>% mutate(id = 1:n()) 
x <- x %>% mutate(id = 1:n()) 
join_overlap_inner(x,y)
#> IRanges object with 4 ranges and 2 metadata columns:
#>           start       end     width |      id.x      id.y
#>       <integer> <integer> <integer> | <integer> <integer>
#>   [1]        11        21        11 |         1         1
#>   [2]        11        21        11 |         1         2
#>   [3]       101       201       101 |         2         3
#>   [4]       101       201       101 |         2         4
join_overlap_left(y,x, suffix = c(".left", ".right"))
#> IRanges object with 5 ranges and 2 metadata columns:
#>           start       end     width |   id.left  id.right
#>       <integer> <integer> <integer> | <integer> <integer>
#>   [1]        10        19        10 |         1         1
#>   [2]        20        21         2 |         2         1
#>   [3]        50       105        56 |         3         2
#>   [4]       100       202       103 |         4         2
#>   [5]         1         5         5 |         5      <NA>
```
