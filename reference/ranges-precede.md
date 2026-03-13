# Find preceding Ranges

Find preceding Ranges

## Usage

``` r
join_precede(x, y, suffix = c(".x", ".y"))

join_precede_right(x, y, suffix = c(".x", ".y"))

join_precede_downstream(x, y, suffix = c(".x", ".y"))
```

## Arguments

- x, y:

  Ranges objects, which ranges in x precede those in y.

- suffix:

  A character vector of length two used to identify metadata columns
  coming from x and y.

## Value

A Ranges object corresponding to the ranges in `y` that are preceded by
the ranges in `x`, all metadata is copied over from the right-hand side
ranges `y`.

## Details

By default `join_precede` will return the ranges in x that come before
the ranges in y and ignore any strand information. The function
`join_precede_right` will find all ranges in y that are on the
right-hand side of the ranges in x ignoring any strand information.
Finally, `join_precede_downstream` will find all ranges in y that are
that are downstream of the ranges in x. On the positive strand this will
result in ranges in y that are right of those in x and on the negative
strand it will result in ranges in y that are left of those in x.

## Examples

``` r
subject <- data.frame(start = c(5,10, 15,20), width = 5, gc = runif(4)) %>%
             as_iranges()
query <- data.frame(start = 2:6, width = 3:7, label = letters[1:5]) %>%
             as_iranges()

join_precede(query, subject)
#> IRanges object with 5 ranges and 2 metadata columns:
#>           start       end     width |       label        gc
#>       <integer> <integer> <integer> | <character> <numeric>
#>   [1]         2         4         3 |           a 0.0433559
#>   [2]         3         6         4 |           b 0.0195993
#>   [3]         4         8         5 |           c 0.0195993
#>   [4]         5        10         6 |           d 0.1990725
#>   [5]         6        12         7 |           e 0.1990725

query  <- data.frame(seqnames = "chr1",
               start = c(11,101),
               end = c(21, 200),
               name = c("a1", "a2"),
               strand = c("+", "-"),
               score = c(1,2)) %>%
           as_granges()
subject <- data.frame(seqnames = "chr1",
                      strand = c("+", "-", "+", "-"),
                      start = c(21,91,101,201),
                      end = c(30,101,110,210),
                      name = paste0("b", 1:4),
                      score = 1:4) %>%
                   as_granges()

join_precede(query, subject)
#> GRanges object with 2 ranges and 4 metadata columns:
#>       seqnames    ranges strand |      name.x   score.x      name.y   score.y
#>          <Rle> <IRanges>  <Rle> | <character> <numeric> <character> <integer>
#>   [1]     chr1     11-21      + |          a1         1          b2         2
#>   [2]     chr1   101-200      - |          a2         2          b4         4
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
join_precede_right(query, subject)
#> GRanges object with 2 ranges and 4 metadata columns:
#>       seqnames    ranges strand |      name.x   score.x      name.y   score.y
#>          <Rle> <IRanges>  <Rle> | <character> <numeric> <character> <integer>
#>   [1]     chr1     11-21      + |          a1         1          b2         2
#>   [2]     chr1   101-200      - |          a2         2          b4         4
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
join_precede_downstream(query, subject)
#> GRanges object with 1 range and 4 metadata columns:
#>       seqnames    ranges strand |      name.x   score.x      name.y   score.y
#>          <Rle> <IRanges>  <Rle> | <character> <numeric> <character> <integer>
#>   [1]     chr1     11-21      + |          a1         1          b3         3
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
