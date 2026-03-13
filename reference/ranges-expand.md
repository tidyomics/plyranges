# Expand list-columns in a Ranges object

Expand list-columns in a Ranges object

## Usage

``` r
expand_ranges(
  data,
  ...,
  .drop = FALSE,
  .id = NULL,
  .keep_empty = FALSE,
  .recursive = FALSE
)
```

## Arguments

- data:

  A Ranges object

- ...:

  list-column names to expand then unlist

- .drop:

  Should additional list columns be dropped (default = FALSE)? By
  default `expand_ranges()` will keep other list columns even if they
  are nested.

- .id:

  A character vector of length equal to number of list columns. If
  supplied will create new column(s) with name `.id` identifying the
  index of the list column (default = NULL).

- .keep_empty:

  If a list-like column contains empty elements, should those elements
  be kept? (default = FALSE)

- .recursive:

  If there are multiple list-columns, should the columns be treated as
  parallel? If FALSE each column will be unnested recursively, otherwise
  they are treated as parallel, that is each list column has identical
  lengths. (deafualt = FALSE)

## Value

a GRanges object with expanded list columns

## Examples

``` r
grng <- as_granges(data.frame(seqnames = "chr1", start = 20:23, width = 1000))
grng <- mutate(grng, 
               exon_id = IntegerList(a = 1, b = c(4,5), c = 3, d = c(2,5))
               )
expand_ranges(grng)
#> GRanges object with 6 ranges and 1 metadata column:
#>       seqnames    ranges strand |   exon_id
#>          <Rle> <IRanges>  <Rle> | <integer>
#>   [1]     chr1   20-1019      * |         1
#>   [2]     chr1   21-1020      * |         4
#>   [3]     chr1   21-1020      * |         5
#>   [4]     chr1   22-1021      * |         3
#>   [5]     chr1   23-1022      * |         2
#>   [6]     chr1   23-1022      * |         5
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
expand_ranges(grng, .id = "name")
#> GRanges object with 6 ranges and 2 metadata columns:
#>       seqnames    ranges strand |   exon_id        name
#>          <Rle> <IRanges>  <Rle> | <integer> <character>
#>   [1]     chr1   20-1019      * |         1           a
#>   [2]     chr1   21-1020      * |         4           b
#>   [3]     chr1   21-1020      * |         5           b
#>   [4]     chr1   22-1021      * |         3           c
#>   [5]     chr1   23-1022      * |         2           d
#>   [6]     chr1   23-1022      * |         5           d
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

# empty list elements are not preserved by default
grng <- mutate(grng, 
               exon_id = IntegerList(a = NULL, b = c(4,5), c= 3, d = c(2,5))
               )
expand_ranges(grng)
#> GRanges object with 5 ranges and 1 metadata column:
#>       seqnames    ranges strand |   exon_id
#>          <Rle> <IRanges>  <Rle> | <integer>
#>   [1]     chr1   21-1020      * |         4
#>   [2]     chr1   21-1020      * |         5
#>   [3]     chr1   22-1021      * |         3
#>   [4]     chr1   23-1022      * |         2
#>   [5]     chr1   23-1022      * |         5
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
expand_ranges(grng, .keep_empty = TRUE)
#> GRanges object with 6 ranges and 1 metadata column:
#>       seqnames    ranges strand |   exon_id
#>          <Rle> <IRanges>  <Rle> | <integer>
#>   [1]     chr1   20-1019      * |      <NA>
#>   [2]     chr1   21-1020      * |         4
#>   [3]     chr1   21-1020      * |         5
#>   [4]     chr1   22-1021      * |         3
#>   [5]     chr1   23-1022      * |         2
#>   [6]     chr1   23-1022      * |         5
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
expand_ranges(grng, .id = "name", .keep_empty = TRUE)
#> GRanges object with 6 ranges and 2 metadata columns:
#>       seqnames    ranges strand |   exon_id        name
#>          <Rle> <IRanges>  <Rle> | <integer> <character>
#>   [1]     chr1   20-1019      * |      <NA>           a
#>   [2]     chr1   21-1020      * |         4           b
#>   [3]     chr1   21-1020      * |         5           b
#>   [4]     chr1   22-1021      * |         3           c
#>   [5]     chr1   23-1022      * |         2           d
#>   [6]     chr1   23-1022      * |         5           d
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
