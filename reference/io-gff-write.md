# Write a GFF(123) file

This is a lightweight wrapper to the export family of functions defined
in rtracklayer.

## Usage

``` r
write_gff(x, file, index = FALSE)

write_gff1(x, file, index = FALSE)

write_gff2(x, file, index = FALSE)

write_gff3(x, file, index = FALSE)
```

## Arguments

- x:

  A GRanges object

- file:

  Path or connection to write to

- index:

  If TRUE the output file will be compressed and indexed using bgzf and
  tabix.

## Value

The write function returns a GFFFile object invisibly

## See also

`rtracklayer::`[`GFFFile()`](https://rdrr.io/pkg/rtracklayer/man/GFFFile-class.html)

## Examples

``` r
if (FALSE) { # \dontrun{
 test_path <- system.file("tests", package = "rtracklayer")
 test_gff3 <- file.path(test_path, "genes.gff3")
 gr <- read_gff3(test_gff3)
 out_gff3 <- file.path(tempdir(), "test.gff3")
 write_gff3(gr, out_gff3)
 read_gff3(out_gff3)
} # }
```
