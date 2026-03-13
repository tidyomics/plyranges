# Pair together two ranges objects

Pair together two ranges objects

## Usage

``` r
pair_overlaps(x, y, maxgap, minoverlap, suffix)

pair_nearest(x, y, suffix)

pair_precede(x, y, suffix)

pair_follow(x, y, suffix)
```

## Arguments

- x, y:

  Ranges objects to pair together.

- maxgap, minoverlap:

  The maximimum gap between intervals as an integer greater than or
  equal to negative one. The minimum amount of overlap between intervals
  as an integer greater than zero, accounting for the maximum gap.

- suffix:

  A character vector of length two used to identify metadata columns
  coming from x and y.

## Value

a DataFrame with two ranges columns and the corresponding metadata
columns.

## Details

These functions return a DataFrame object, and is one way of
representing paired alignments with plyranges.

## See also

\[join_nearest()\]\[join_overlap_inner()\]\[join_precede()\]\[join_follow()\]

## Examples

``` r
query <- data.frame(start = c(5,10, 15,20), width = 5, gc = runif(4)) %>%
             as_iranges()
subject <- data.frame(start = 2:6, width = 3:7, label = letters[1:5]) %>%
             as_iranges()

pair_overlaps(query, subject)
#> DataFrame with 6 rows and 4 columns
#>    ranges.x  ranges.y        gc       label
#>   <IRanges> <IRanges> <numeric> <character>
#> 1       5-9       3-6  0.908478           b
#> 2       5-9       4-8  0.908478           c
#> 3       5-9      5-10  0.908478           d
#> 4       5-9      6-12  0.908478           e
#> 5     10-14      5-10  0.122533           d
#> 6     10-14      6-12  0.122533           e
pair_overlaps(query, subject, minoverlap = 5)
#> DataFrame with 1 row and 4 columns
#>    ranges.x  ranges.y        gc       label
#>   <IRanges> <IRanges> <numeric> <character>
#> 1       5-9      5-10  0.908478           d
pair_nearest(query, subject)
#> DataFrame with 4 rows and 4 columns
#>    ranges.x  ranges.y        gc       label
#>   <IRanges> <IRanges> <numeric> <character>
#> 1       5-9      6-12  0.908478           e
#> 2     10-14      6-12  0.122533           e
#> 3     15-19      6-12  0.728863           e
#> 4     20-24      6-12  0.950377           e


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
pair_overlaps(query, subject, suffix = c(".query", ".subject"))
#> DataFrame with 3 rows and 6 columns
#>    granges.query granges.subject  name.query score.query name.subject
#>        <GRanges>       <GRanges> <character>   <numeric>  <character>
#> 1   chr1:11-21:+    chr1:21-30:+          a1           1           b1
#> 2 chr1:101-200:-   chr1:91-101:-          a2           2           b2
#> 3 chr1:101-200:-  chr1:101-110:+          a2           2           b3
#>   score.subject
#>       <integer>
#> 1             1
#> 2             2
#> 3             3
pair_follow(query, subject, suffix = c(".query", ".subject"))
#> DataFrame with 1 row and 6 columns
#>    granges.query granges.subject  name.query score.query name.subject
#>        <GRanges>       <GRanges> <character>   <numeric>  <character>
#> 1 chr1:101-200:-    chr1:21-30:+          a2           2           b1
#>   score.subject
#>       <integer>
#> 1             1
pair_precede(query, subject, suffix = c(".query", ".subject"))
#> DataFrame with 2 rows and 6 columns
#>    granges.query granges.subject  name.query score.query name.subject
#>        <GRanges>       <GRanges> <character>   <numeric>  <character>
#> 1   chr1:11-21:+   chr1:91-101:-          a1           1           b2
#> 2 chr1:101-200:-  chr1:201-210:-          a2           2           b4
#>   score.subject
#>       <integer>
#> 1             2
#> 2             4
pair_precede(query, subject, suffix = c(".query", ".subject"))
#> DataFrame with 2 rows and 6 columns
#>    granges.query granges.subject  name.query score.query name.subject
#>        <GRanges>       <GRanges> <character>   <numeric>  <character>
#> 1   chr1:11-21:+   chr1:91-101:-          a1           1           b2
#> 2 chr1:101-200:-  chr1:201-210:-          a2           2           b4
#>   score.subject
#>       <integer>
#> 1             2
#> 2             4
```
