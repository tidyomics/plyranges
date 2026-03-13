# Anchored Ranges objects

The `GRangesAnchored` class and the `IRangesAnchored` class allow
components of a `GRanges` or `IRanges` (start, end, center) to be held
fixed.

## Usage

``` r
anchor(x)

unanchor(x)

anchor_start(x)

anchor_end(x)

anchor_center(x)

anchor_centre(x)

anchor_3p(x)

anchor_5p(x)
```

## Arguments

- x:

  a Ranges object

## Value

a RangesAnchored object which has the same appearance as a regular
Ranges object but with an additional slot displaying an anchor.

## Details

Anchoring will fix a Ranges start, end, or center positions, so these
positions will remain the same when performing arithimetic. For
`GRanges` objects, the function (`anchor_3p()`) will fix the start for
the negative strand, while `anchor_5p()` will fix the end for the
positive strand. Anchoring modifies how arithmetic is performed, for
example modifying the width of a range with
[`set_width()`](https://tidyomics.github.io/plyranges/reference/ranges-setters.md)
or stretching a range with
[`stretch()`](https://tidyomics.github.io/plyranges/reference/stretch.md).
To remove anchoring use `unanchor()`.

## Constructors

Depending on how you want to fix the components of a Ranges, there are
five ways to construct a RangesAnchored class. Here `x` is either an
`IRanges` or `GRanges` object.

- `anchor_start(x)`:

  Fix the start coordinates

- `anchor_end(x)`:

  Fix the end coordinates

- `anchor_center(x)`:

  Fix the center coordinates

- `anchor_3p(x)`:

  On the negative strand fix the start coordinates, and for positive or
  unstranded ranges fix the end coordinates.

- `anchor_5p(x)`:

  On the positive or unstranded ranges fix the start coordinates,
  coordinates and for negative stranded ranges fix the end coordinates.

## Accessors

To see what has been anchored use the function `anchor`. This will
return a character vector containing a valid anchor. It will be set to
one of `c("start", "end", "center")` for an `IRanges` object or one of
`c("start", "end", "center", "3p", "5p")` for a `GRanges` object.

## See also

[mutate.Ranges](https://tidyomics.github.io/plyranges/reference/mutate-ranges.md),
[stretch](https://tidyomics.github.io/plyranges/reference/stretch.md)

## Examples

``` r
df <- data.frame(start = 1:10, width = 5)
rng <- as_iranges(df)
rng_by_start <- anchor_start(rng)
rng_by_start
#> IRanges object with 10 ranges and 0 metadata columns:
#> Anchored by: start
#>            start       end     width
#>        <integer> <integer> <integer>
#>    [1]         1         5         5
#>    [2]         2         6         5
#>    [3]         3         7         5
#>    [4]         4         8         5
#>    [5]         5         9         5
#>    [6]         6        10         5
#>    [7]         7        11         5
#>    [8]         8        12         5
#>    [9]         9        13         5
#>   [10]        10        14         5
anchor(rng_by_start)
#> [1] "start"
mutate(rng_by_start, width = 3L)
#> IRanges object with 10 ranges and 0 metadata columns:
#>            start       end     width
#>        <integer> <integer> <integer>
#>    [1]         1         3         3
#>    [2]         2         4         3
#>    [3]         3         5         3
#>    [4]         4         6         3
#>    [5]         5         7         3
#>    [6]         6         8         3
#>    [7]         7         9         3
#>    [8]         8        10         3
#>    [9]         9        11         3
#>   [10]        10        12         3
grng <- as_granges(df,
                   seqnames = "chr1",
                   strand = c(rep("-", 5), rep("+", 5)))
rng_by_5p <- anchor_5p(grng)
rng_by_5p
#> GRanges object with 10 ranges and 0 metadata columns:
#> Anchored by: 5p
#>        seqnames    ranges strand
#>           <Rle> <IRanges>  <Rle>
#>    [1]     chr1       1-5      -
#>    [2]     chr1       2-6      -
#>    [3]     chr1       3-7      -
#>    [4]     chr1       4-8      -
#>    [5]     chr1       5-9      -
#>    [6]     chr1      6-10      +
#>    [7]     chr1      7-11      +
#>    [8]     chr1      8-12      +
#>    [9]     chr1      9-13      +
#>   [10]     chr1     10-14      +
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
mutate(rng_by_5p, width = 3L)
#> GRanges object with 10 ranges and 0 metadata columns:
#>        seqnames    ranges strand
#>           <Rle> <IRanges>  <Rle>
#>    [1]     chr1       3-5      -
#>    [2]     chr1       4-6      -
#>    [3]     chr1       5-7      -
#>    [4]     chr1       6-8      -
#>    [5]     chr1       7-9      -
#>    [6]     chr1       6-8      +
#>    [7]     chr1       7-9      +
#>    [8]     chr1      8-10      +
#>    [9]     chr1      9-11      +
#>   [10]     chr1     10-12      +
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
