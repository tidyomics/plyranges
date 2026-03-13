# Find following Ranges

Find following Ranges

## Usage

``` r
join_follow(x, y, suffix = c(".x", ".y"))

join_follow_left(x, y, suffix = c(".x", ".y"))

join_follow_upstream(x, y, suffix = c(".x", ".y"))
```

## Arguments

- x, y:

  Ranges objects, which ranges in x follow those in y.

- suffix:

  A character vector of length two used to identify metadata columns
  coming from x and y.

## Value

A Ranges object corresponding to the ranges in
``` x`` that are followed by the ranges in  ```y`, all metadata is copied over from the right-hand side ranges `y\`.

## Details

By default `join_follow` will find abritrary ranges in y that are
followed by ranges in x and ignore any strand information. On the other
hand `join_follow_left` will find all ranges in y that are on the
left-hand side of the ranges in x ignoring any strand information.
Finally, `join_follow_upstream` will find all ranges in x that are that
are upstream of the ranges in y. On the positive strand this will result
in ranges in y that are left of those in x and on the negative strand it
will result in ranges in y that are right of those in x.

## Examples

``` r
query <- data.frame(start = c(5,10, 15,20), width = 5, gc = runif(4)) %>%
             as_iranges()
subject <- data.frame(start = 2:6, width = 3:7, label = letters[1:5]) %>%
             as_iranges()

join_follow(query, subject)
#> IRanges object with 4 ranges and 2 metadata columns:
#>           start       end     width |        gc       label
#>       <integer> <integer> <integer> | <numeric> <character>
#>   [1]         5         9         5 |  0.957558           a
#>   [2]        10        14         5 |  0.551651           c
#>   [3]        15        19         5 |  0.101954           e
#>   [4]        20        24         5 |  0.237915           e

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

join_follow(query, subject)
#> GRanges object with 3 ranges and 4 metadata columns:
#>       seqnames    ranges strand |      name.x   score.x      name.y   score.y
#>          <Rle> <IRanges>  <Rle> | <character> <integer> <character> <numeric>
#>   [1]     chr1    91-101      - |          b2         2          a1         1
#>   [2]     chr1   101-110      + |          b3         3          a1         1
#>   [3]     chr1   201-210      - |          b4         4          a2         2
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
join_follow_left(query, subject)
#> GRanges object with 3 ranges and 4 metadata columns:
#>       seqnames    ranges strand |      name.x   score.x      name.y   score.y
#>          <Rle> <IRanges>  <Rle> | <character> <integer> <character> <numeric>
#>   [1]     chr1    91-101      - |          b2         2          a1         1
#>   [2]     chr1   101-110      + |          b3         3          a1         1
#>   [3]     chr1   201-210      - |          b4         4          a2         2
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
join_follow_upstream(query, subject)
#> GRanges object with 1 range and 4 metadata columns:
#>       seqnames    ranges strand |      name.x   score.x      name.y   score.y
#>          <Rle> <IRanges>  <Rle> | <character> <integer> <character> <numeric>
#>   [1]     chr1   101-110      + |          b3         3          a1         1
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
