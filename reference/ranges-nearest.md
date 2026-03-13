# Find nearest neighbours between two Ranges objects

Find nearest neighbours between two Ranges objects

## Usage

``` r
join_nearest(x, y, suffix = c(".x", ".y"), distance = FALSE)

join_nearest_left(x, y, suffix = c(".x", ".y"), distance = FALSE)

join_nearest_right(x, y, suffix = c(".x", ".y"), distance = FALSE)

join_nearest_upstream(x, y, suffix = c(".x", ".y"), distance = FALSE)

join_nearest_downstream(x, y, suffix = c(".x", ".y"), distance = FALSE)
```

## Arguments

- x, y:

  Ranges objects, add the nearest neighbours of ranges in x to those in
  y.

- suffix:

  A character vector of length two used to identify metadata columns

- distance:

  logical vector whether to add a column named "distance" containing the
  distance to the nearest region. If set to a character vector of length
  1, will use that as distance column name.

## Value

A Ranges object corresponding to the nearest ranges, all metadata is
copied over from the right-hand side ranges `y`.

## Details

By default `join_nearest` will find arbitrary nearest neighbours in
either direction and ignore any strand information. The
`join_nearest_left` and `join_nearest_right` methods will find arbitrary
nearest neighbour ranges on x that are left/right of those on y and
ignore any strand information.

The `join_nearest_upstream` method will find arbitrary nearest neighbour
ranges on x that are upstream of those on y. This takes into account
strandedness of the ranges. On the positive strand nearest upstream will
be on the left and on the negative strand nearest upstream will be on
the right.

The `join_nearest_downstream` method will find arbitrary nearest
neighbour ranges on x that are upstream of those on y. This takes into
account strandedness of the ranges.On the positive strand nearest
downstream will be on the right and on the negative strand nearest
upstream will be on the left.

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

join_nearest(query, subject)
#> IRanges object with 4 ranges and 2 metadata columns:
#>           start       end     width |        gc       label
#>       <integer> <integer> <integer> | <numeric> <character>
#>   [1]         5         9         5 |  0.860098           e
#>   [2]        10        14         5 |  0.738261           e
#>   [3]        15        19         5 |  0.497312           e
#>   [4]        20        24         5 |  0.579921           f
join_nearest_left(query, subject)
#> IRanges object with 3 ranges and 2 metadata columns:
#>           start       end     width |        gc       label
#>       <integer> <integer> <integer> | <numeric> <character>
#>   [1]         5         9         5 |  0.860098           a
#>   [2]        10        14         5 |  0.738261           d
#>   [3]        15        19         5 |  0.497312           e
join_nearest_right(query, subject)
#> IRanges object with 1 range and 2 metadata columns:
#>           start       end     width |        gc       label
#>       <integer> <integer> <integer> | <numeric> <character>
#>   [1]        20        24         5 |  0.579921           f

subject  <- data.frame(seqnames = "chr1",
               start = c(11,101),
               end = c(21, 200),
               name = c("a1", "a2"),
               strand = c("+", "-"),
               score = c(1,2)) %>%
           as_granges()
query <- data.frame(seqnames = "chr1",
                      strand = c("+", "-", "+", "-"),
                      start = c(21,91,101,201),
                      end = c(30,101,110,210),
                      name = paste0("b", 1:4),
                      score = 1:4) %>%
                   as_granges()
join_nearest_upstream(query, subject)
#> GRanges object with 3 ranges and 4 metadata columns:
#>       seqnames    ranges strand |      name.x   score.x      name.y   score.y
#>          <Rle> <IRanges>  <Rle> | <character> <integer> <character> <numeric>
#>   [1]     chr1     21-30      + |          b1         1          a1         1
#>   [2]     chr1    91-101      - |          b2         2          a2         2
#>   [3]     chr1   101-110      + |          b3         3          a1         1
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
join_nearest_downstream(query, subject)
#> GRanges object with 1 range and 4 metadata columns:
#>       seqnames    ranges strand |      name.x   score.x      name.y   score.y
#>          <Rle> <IRanges>  <Rle> | <character> <integer> <character> <numeric>
#>   [1]     chr1   201-210      - |          b4         4          a2         2
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
