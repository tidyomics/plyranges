# Write a BigWig file

This is a lightweight wrapper to the export family of functions defined
in rtracklayer.

## Usage

``` r
write_bigwig(x, file)
```

## Arguments

- x:

  A GRanges object

- file:

  File name, URL or connection specifying a file to write x to.
  Compressed files with extensions such as '.gz' are handled
  automatically.

## Value

The write functions return a BigWigFile invisibly

## See also

`rtracklayer::`[`BigWigFile()`](https://rdrr.io/pkg/rtracklayer/man/BigWigFile.html)

## Examples

``` r
if (FALSE) { # \dontrun{
 if (.Platform$OS.type != "windows") {
  test_path <- system.file("tests", package = "rtracklayer")
  bw_file <- file.path(test_path, "test.bw")
  gr <- read_bigwig(bw_file)
  gr
  bw_out <- file.path(tempdir(), "test_out.bw")
  write_bigwig(gr ,bw_out)
  read_bigwig(bw_out)
 }
} # }
```
