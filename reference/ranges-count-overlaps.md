# Count the number of overlaps between two Ranges objects

Count the number of overlaps between two Ranges objects

## Usage

``` r
count_overlaps(x, y, maxgap, minoverlap)

# S3 method for class 'IntegerRanges'
count_overlaps(x, y, maxgap = -1L, minoverlap = 0L)

# S3 method for class 'GenomicRanges'
count_overlaps(x, y, maxgap = -1L, minoverlap = 0L)

count_overlaps_within(x, y, maxgap, minoverlap)

# S3 method for class 'IntegerRanges'
count_overlaps_within(x, y, maxgap = 0L, minoverlap = 1L)

# S3 method for class 'GenomicRanges'
count_overlaps_within(x, y, maxgap = 0L, minoverlap = 1L)

count_overlaps_directed(x, y, maxgap, minoverlap)

# S3 method for class 'GenomicRanges'
count_overlaps_directed(x, y, maxgap = -1L, minoverlap = 0L)

count_overlaps_within_directed(x, y, maxgap, minoverlap)

# S3 method for class 'GenomicRanges'
count_overlaps_within_directed(x, y, maxgap = -1L, minoverlap = 0L)
```

## Arguments

- x, y:

  Objects representing ranges

- maxgap, minoverlap:

  The maximimum gap between intervals as an integer greater than or
  equal to zero. The minimum amount of overlap between intervals as an
  integer greater than zero, accounting for the maximum gap.

## Value

An integer vector of same length as x.

## Examples

``` r
query <- data.frame(start = c(5,10, 15,20), width = 5, gc = runif(4)) %>%
             as_iranges()
subject <- data.frame(start = 2:6, width = 3:7, label = letters[1:5]) %>%
             as_iranges()
query %>% mutate(n_olap = count_overlaps(., subject),
                 n_olap_within = count_overlaps_within(., subject))
#> Error: unable to find an inherited method for function ‘countOverlaps’ for signature ‘query = "IRanges", subject = "standardGeneric"’
```
