# Construct annotation information

To construct annotations by supplying annotation information use
`genome_info`. To add annotations to an existing Ranges object use
`set_genome_info`. To retrieve an annotation as a Ranges object use
`get_genome_info`.

## Usage

``` r
genome_info(
  genome = NULL,
  seqnames = NULL,
  seqlengths = NULL,
  is_circular = NULL
)

set_genome_info(
  .data,
  genome = NULL,
  seqnames = NULL,
  seqlengths = NULL,
  is_circular = NULL
)

get_genome_info(.data)
```

## Arguments

- genome:

  A character vector of length one indicating the genome build.

- seqnames:

  A character vector containing the name of sequences.

- seqlengths:

  An optional integer vector containg the lengths of sequences.

- is_circular:

  An optional logical vector indicating whether a sequence is ciruclar.

- .data:

  A Ranges object to annotate or retrieve an annotation for.

## Value

a GRanges object containing annotations. To retrieve the annotations as
a Ranges object use `get_genome_info`.

## See also

[`Seqinfo::Seqinfo()`](https://rdrr.io/pkg/Seqinfo/man/Seqinfo-class.html)

## Examples

``` r
x <- genome_info(genome = "toy",
                 seqnames = letters[1:4],
                 seqlengths = c(100, 300, 15, 600),
                 is_circular = c(NA, FALSE, FALSE, TRUE))
x
#> GRanges object with 4 ranges and 0 metadata columns:
#>     seqnames    ranges strand
#>        <Rle> <IRanges>  <Rle>
#>   a        a     1-100      *
#>   b        b     1-300      *
#>   c        c      1-15      *
#>   d        d     1-600      *
#>   -------
#>   seqinfo: 4 sequences (1 circular) from toy genome

rng <- as_granges(data.frame(seqnames = "a", start = 30:50, width = 10))
rng
#> GRanges object with 21 ranges and 0 metadata columns:
#>        seqnames    ranges strand
#>           <Rle> <IRanges>  <Rle>
#>    [1]        a     30-39      *
#>    [2]        a     31-40      *
#>    [3]        a     32-41      *
#>    [4]        a     33-42      *
#>    [5]        a     34-43      *
#>    ...      ...       ...    ...
#>   [17]        a     46-55      *
#>   [18]        a     47-56      *
#>   [19]        a     48-57      *
#>   [20]        a     49-58      *
#>   [21]        a     50-59      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
rng <- set_genome_info(rng,
                       genome = "toy",
                       seqnames = letters[1:4],
                       seqlengths = c(100, 300, 15, 600),
                       is_circular = c(NA, FALSE, FALSE, TRUE))
get_genome_info(rng)
#> GRanges object with 4 ranges and 0 metadata columns:
#>     seqnames    ranges strand
#>        <Rle> <IRanges>  <Rle>
#>   a        a     1-100      *
#>   b        b     1-300      *
#>   c        c      1-15      *
#>   d        d     1-600      *
#>   -------
#>   seqinfo: 4 sequences (1 circular) from toy genome

if (FALSE) { # \dontrun{
if (interactive()) {
 # requires internet connection
 genome_info(genome = "hg38")
}
} # }
```
