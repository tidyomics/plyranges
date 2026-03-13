# Read a WIG file

This is a lightweight wrapper to the import family of functions defined
in rtracklayer.

## Usage

``` r
read_wig(file, genome_info = NULL, overlap_ranges = NULL)
```

## Arguments

- file:

  A path to a file or a connection.

- genome_info:

  An optional character string or a Ranges object that contains
  information about the genome build. For example the USSC identifier
  "hg19" will add build information to the returned GRanges.

- overlap_ranges:

  An optional Ranges object. Only the intervals in the file that overlap
  the Ranges will be returned.

## Value

A GRanges object

A GRanges object

## See also

`rtracklayer::`[`WIGFile()`](https://rdrr.io/pkg/rtracklayer/man/WIGFile-class.html)

## Examples

``` r
test_path <- system.file("tests", package = "rtracklayer")
test_wig <- file.path(test_path, "step.wig")
gr <- read_wig(test_wig)
gr
#> GRanges object with 19 ranges and 1 metadata column:
#>        seqnames            ranges strand |     score
#>           <Rle>         <IRanges>  <Rle> | <numeric>
#>    [1]    chr19          59104701      * |      10.0
#>    [2]    chr19          59104901      * |      12.5
#>    [3]    chr19          59105401      * |      15.0
#>    [4]    chr19          59105601      * |      17.5
#>    [5]    chr19          59105901      * |      20.0
#>    ...      ...               ...    ... .       ...
#>   [15]    chr18 59109521-59109720      * |       500
#>   [16]    chr18 59109821-59110020      * |       400
#>   [17]    chr18 59110121-59110320      * |       300
#>   [18]    chr18 59110421-59110620      * |       200
#>   [19]    chr18 59110721-59110920      * |       100
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
gr <- read_wig(test_wig, genome_info = "hg19")
```
