# Find overlap between two Ranges

Find overlap between two Ranges

## Usage

``` r
find_overlaps(x, y, maxgap, minoverlap, suffix = c(".x", ".y"))

# S3 method for class 'IntegerRanges'
find_overlaps(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y"))

# S3 method for class 'GenomicRanges'
find_overlaps(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y"))

find_overlaps_within(x, y, maxgap, minoverlap, suffix = c(".x", ".y"))

# S3 method for class 'IntegerRanges'
find_overlaps_within(
  x,
  y,
  maxgap = -1L,
  minoverlap = 0L,
  suffix = c(".x", ".y")
)

# S3 method for class 'GenomicRanges'
find_overlaps_within(
  x,
  y,
  maxgap = -1L,
  minoverlap = 0L,
  suffix = c(".x", ".y")
)

find_overlaps_directed(x, y, maxgap, minoverlap, suffix = c(".x", ".y"))

# S3 method for class 'GenomicRanges'
find_overlaps_directed(
  x,
  y,
  maxgap = -1L,
  minoverlap = 0L,
  suffix = c(".x", ".y")
)

find_overlaps_within_directed(x, y, maxgap, minoverlap, suffix = c(".x", ".y"))

# S3 method for class 'GenomicRanges'
find_overlaps_within_directed(x, y, maxgap, minoverlap, suffix = c(".x", ".y"))

group_by_overlaps(x, y, maxgap, minoverlap)

# S3 method for class 'IntegerRanges'
group_by_overlaps(x, y, maxgap = -1L, minoverlap = 0L)

# S3 method for class 'GenomicRanges'
group_by_overlaps(x, y, maxgap = -1L, minoverlap = 0L)
```

## Arguments

- x, y:

  Objects representing ranges

- maxgap, minoverlap:

  The maximimum gap between intervals as an integer greater than or
  equal to negative one. The minimum amount of overlap between intervals
  as an integer greater than zero, accounting for the maximum gap.

- suffix:

  A character vector of length two used to identify metadata columns
  coming from x and y.

## Value

A Ranges object with rows corresponding to the ranges in x that overlap
y. In the case of `group_by_overlaps()`, returns a GroupedRanges object,
grouped by the number of overlaps of ranges in x that overlap y (stored
in a column called query).

## Details

`find_overlaps()` will search for any overlaps between ranges x and y
and return a Ranges object of length equal to the number of times x
overlaps y. This Ranges object will have additional metadata columns
corresponding to the metadata columns in y. `find_overlaps_within()` is
the same but will only search for overlaps within y. For GRanges objects
strand is ignored, unless `find_overlaps_directed()` is used. If the
Ranges objects have no metadata, one could use `group_by_overlaps()` to
be able to identify the index of the input Range x that overlaps a Range
in y. Alternatively,
[`pair_overlaps()`](https://tidyomics.github.io/plyranges/reference/ranges-pairs.md)
could be used to place the x ranges next to the range in y they overlap.

## See also

`IRanges::`[`findOverlaps()`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html),
`GenomicRanges::`[`findOverlaps()`](https://rdrr.io/pkg/GenomicRanges/man/findOverlaps-methods.html)

## Examples

``` r
query <- data.frame(start = c(5,10, 15,20), width = 5, gc = runif(4)) %>%
             as_iranges()
subject <- data.frame(start = 2:6, width = 3:7, label = letters[1:5]) %>%
             as_iranges()

find_overlaps(query, subject)
#> IRanges object with 6 ranges and 2 metadata columns:
#>           start       end     width |        gc       label
#>       <integer> <integer> <integer> | <numeric> <character>
#>   [1]         5         9         5 |  0.629995           b
#>   [2]         5         9         5 |  0.629995           c
#>   [3]         5         9         5 |  0.629995           d
#>   [4]         5         9         5 |  0.629995           e
#>   [5]        10        14         5 |  0.673249           d
#>   [6]        10        14         5 |  0.673249           e
find_overlaps(query, subject, minoverlap = 5)
#> IRanges object with 1 range and 2 metadata columns:
#>           start       end     width |        gc       label
#>       <integer> <integer> <integer> | <numeric> <character>
#>   [1]         5         9         5 |  0.629995           d
find_overlaps_within(query, subject) # same result as minoverlap
#> IRanges object with 1 range and 2 metadata columns:
#>           start       end     width |        gc       label
#>       <integer> <integer> <integer> | <numeric> <character>
#>   [1]         5         9         5 |  0.629995           d
find_overlaps(query, subject, maxgap = 1)
#> IRanges object with 8 ranges and 2 metadata columns:
#>           start       end     width |        gc       label
#>       <integer> <integer> <integer> | <numeric> <character>
#>   [1]         5         9         5 |  0.629995           a
#>   [2]         5         9         5 |  0.629995           b
#>   [3]         5         9         5 |  0.629995           c
#>   [4]         5         9         5 |  0.629995           d
#>   [5]         5         9         5 |  0.629995           e
#>   [6]        10        14         5 |  0.673249           c
#>   [7]        10        14         5 |  0.673249           d
#>   [8]        10        14         5 |  0.673249           e

# -- GRanges objects, strand is ignored by default
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

# ignores strandedness
find_overlaps(query, subject, suffix = c(".query", ".subject"))
#> GRanges object with 3 ranges and 4 metadata columns:
#>       seqnames    ranges strand |  name.query score.query name.subject
#>          <Rle> <IRanges>  <Rle> | <character>   <numeric>  <character>
#>   [1]     chr1     11-21      + |          a1           1           b1
#>   [2]     chr1   101-200      - |          a2           2           b2
#>   [3]     chr1   101-200      - |          a2           2           b3
#>       score.subject
#>           <integer>
#>   [1]             1
#>   [2]             2
#>   [3]             3
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
find_overlaps(query, subject, suffix = c(".query", ".subject"), minoverlap = 2)
#> GRanges object with 1 range and 4 metadata columns:
#>       seqnames    ranges strand |  name.query score.query name.subject
#>          <Rle> <IRanges>  <Rle> | <character>   <numeric>  <character>
#>   [1]     chr1   101-200      - |          a2           2           b3
#>       score.subject
#>           <integer>
#>   [1]             3
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
# adding directed prefix includes strand
find_overlaps_directed(query, subject, suffix = c(".query", ".subject"))
#> GRanges object with 2 ranges and 4 metadata columns:
#>       seqnames    ranges strand |  name.query score.query name.subject
#>          <Rle> <IRanges>  <Rle> | <character>   <numeric>  <character>
#>   [1]     chr1     11-21      + |          a1           1           b1
#>   [2]     chr1   101-200      - |          a2           2           b2
#>       score.subject
#>           <integer>
#>   [1]             1
#>   [2]             2
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
