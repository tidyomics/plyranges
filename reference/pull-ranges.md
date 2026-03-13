# Extract a single column from a Ranges object as a vector

Extract a single column from a Ranges object as a vector

## Usage

``` r
# S3 method for class 'Ranges'
pull(.data, var = -1, name = NULL, ...)
```

## Arguments

- .data:

  a `Ranges` object

- var:

  A variable specified as:

  - a literal variable name

  - a positive integer, giving the position counting from the left. In
    this case order is start, end, width, (strand, seqnames), gc and
    score.

  - a negative integer, giving the position counting from the right. The
    default returns the last column (on the assumption that's the column
    you've created most recently). This argument is taken by expression
    and supports quasiquotation (you can unquote column names and column
    locations).

- name:

  An optional parameter that specifies the column to be used as names
  for a named vector. Specified in a similar manner as `var`.

- ...:

  For use by methods.

## See also

[`dplyr::pull()`](https://dplyr.tidyverse.org/reference/pull.html)

## Examples

``` r
df <- data.frame(start = 1:10,
                 width = 5,
                 seqnames = "seq1",
                 strand = sample(c("+", "-", "*"), 10, replace = TRUE),
                 gc = runif(10),
                 score = rpois(10, 2))
rng <- as_granges(df)

# Pull parts of the range
pull(rng, start)
#>  [1]  1  2  3  4  5  6  7  8  9 10
# equivalent to start(rng)

# Pull by column name
pull(rng, gc)
#>  [1] 0.6943507 0.4122370 0.3277259 0.5725648 0.9669991 0.6617790 0.6246977
#>  [8] 0.8566530 0.7747789 0.8340271
pull(rng, score)
#>  [1] 0 2 2 4 5 0 2 3 1 1

# Pull by position (positive from left, negative from right)
pull(rng, 1)    # First metadata column
#>  [1]  1  2  3  4  5  6  7  8  9 10
pull(rng, -1)   # Last metadata column (default)
#>  [1] 0 2 2 4 5 0 2 3 1 1
pull(rng, -2)   # Second to last metadata column
#>  [1] 0.6943507 0.4122370 0.3277259 0.5725648 0.9669991 0.6617790 0.6246977
#>  [8] 0.8566530 0.7747789 0.8340271

# Pull with names
pull(rng, score, name = gc)
#> 0.694350712001324 0.412237035110593 0.327725868672132  0.57256476697512 
#>                 0                 2                 2                 4 
#> 0.966999084455892 0.661779022077098 0.624697716906667 0.856653042603284 
#>                 5                 0                 2                 3 
#>  0.77477888809517 0.834027098724619 
#>                 1                 1 
```
