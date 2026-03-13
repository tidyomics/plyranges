# Read a BigWig file

Read a BigWig file

## Usage

``` r
read_bigwig(file, genome_info = NULL, overlap_ranges = NULL)
```

## Arguments

- file:

  A path to a file or URL.

- genome_info:

  An optional character string or a Ranges object that contains
  information about the genome build. For example the identifier "hg19"
  will add build information to the returned GRanges.

- overlap_ranges:

  An optional Ranges object. Only the intervals in the file that overlap
  the Ranges will be loaded.

## Value

a GRanges object

## See also

`rtracklayer::`[`BigWigFile()`](https://rdrr.io/pkg/rtracklayer/man/BigWigFile.html)

## Examples

``` r
if (.Platform$OS.type != "windows") {
  test_path <- system.file("tests", package = "rtracklayer")
  bw_file <- file.path(test_path, "test.bw")
  gr <- read_bigwig(bw_file)
  gr
}
#> GRanges object with 9 ranges and 1 metadata column:
#>       seqnames    ranges strand |     score
#>          <Rle> <IRanges>  <Rle> | <numeric>
#>   [1]     chr2     1-300      * |     -1.00
#>   [2]     chr2   301-600      * |     -0.75
#>   [3]     chr2   601-900      * |     -0.50
#>   [4]     chr2  901-1200      * |     -0.25
#>   [5]     chr2 1201-1500      * |      0.00
#>   [6]    chr19 1501-1800      * |      0.25
#>   [7]    chr19 1801-2100      * |      0.50
#>   [8]    chr19 2101-2400      * |      0.75
#>   [9]    chr19 2401-2700      * |      1.00
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome
```
