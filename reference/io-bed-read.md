# Read a BED or BEDGraph file

This is a lightweight wrapper to the import family of functions defined
in rtracklayer.

Read common interval based formats as GRanges.

## Usage

``` r
read_bed(file, col_names = NULL, genome_info = NULL, overlap_ranges = NULL)

read_bed_graph(
  file,
  col_names = NULL,
  genome_info = NULL,
  overlap_ranges = NULL
)

read_narrowpeaks(
  file,
  col_names = NULL,
  genome_info = NULL,
  overlap_ranges = NULL
)
```

## Arguments

- file:

  A path to a file or a connection.

- col_names:

  An optional character vector for including additional columns in
  `file` that are not part of the BED/narrowPeaks specification.

- genome_info:

  An optional character string or a Ranges object that contains
  information about the genome build. For example the USSC identifier
  "hg19" will add build information to the returned GRanges.

- overlap_ranges:

  An optional Ranges object. Only the intervals in the file that overlap
  the Ranges will be returned.

## Value

A GRanges object

## Details

This is a lightweight wrapper to the import family of functions defined
in rtracklayer. The `read_narrowpeaks` function parses the ENCODE
narrowPeak BED format (see
<https://genome.ucsc.edu/FAQ/FAQformat.html#format12> for details.). As
such the parser expects four additional columns called (corresponding to
the narrowPeaks spec):

- signalValue

- pValue

- qValue

- peak

## See also

`rtracklayer::`[`BEDFile()`](https://rdrr.io/pkg/rtracklayer/man/BEDFile-class.html)

## Examples

``` r
test_path <- system.file("tests", package = "rtracklayer")
bed_file <- file.path(test_path, "test.bed")
gr <- read_bed(bed_file)
gr
#> GRanges object with 5 ranges and 5 metadata columns:
#>       seqnames              ranges strand |        name     score     itemRgb
#>          <Rle>           <IRanges>  <Rle> | <character> <numeric> <character>
#>   [1]     chr7 127471197-127472363      + |        Pos1         0     #FF0000
#>   [2]     chr7 127472364-127473530      + |        Pos2         2     #FF0000
#>   [3]     chr7 127473531-127474697      - |        Neg1         0     #FF0000
#>   [4]     chr9 127474698-127475864      + |        Pos3         5     #FF0000
#>   [5]     chr9 127475865-127477031      - |        Neg2         5     #0000FF
#>                     thick                  blocks
#>                 <IRanges>           <IRangesList>
#>   [1] 127471197-127472363 1-300,501-700,1068-1167
#>   [2] 127472364-127473530          1-250,668-1167
#>   [3] 127473531-127474697                  1-1167
#>   [4] 127474698-127475864                  1-1167
#>   [5] 127475865-127477031                  1-1167
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
gr <- read_bed(bed_file, genome_info = "hg19")
#> 
#> Attaching package: ‘Biostrings’
#> The following object is masked from ‘package:base’:
#> 
#>     strsplit
gr
#> GRanges object with 5 ranges and 5 metadata columns:
#>       seqnames              ranges strand |        name     score     itemRgb
#>          <Rle>           <IRanges>  <Rle> | <character> <numeric> <character>
#>   [1]     chr7 127471197-127472363      + |        Pos1         0     #FF0000
#>   [2]     chr7 127472364-127473530      + |        Pos2         2     #FF0000
#>   [3]     chr7 127473531-127474697      - |        Neg1         0     #FF0000
#>   [4]     chr9 127474698-127475864      + |        Pos3         5     #FF0000
#>   [5]     chr9 127475865-127477031      - |        Neg2         5     #0000FF
#>                     thick                  blocks
#>                 <IRanges>           <IRangesList>
#>   [1] 127471197-127472363 1-300,501-700,1068-1167
#>   [2] 127472364-127473530          1-250,668-1167
#>   [3] 127473531-127474697                  1-1167
#>   [4] 127474698-127475864                  1-1167
#>   [5] 127475865-127477031                  1-1167
#>   -------
#>   seqinfo: 298 sequences (2 circular) from hg19 genome
olap <-  as_granges(data.frame(seqnames = "chr7", start = 1, end = 127473000))
gr <- read_bed(bed_file,
              overlap_ranges = olap)
# bedGraph
bg_file <- file.path(test_path, "test.bedGraph")
gr <- read_bed_graph(bg_file)
gr
#> GRanges object with 9 ranges and 1 metadata column:
#>       seqnames            ranges strand |     score
#>          <Rle>         <IRanges>  <Rle> | <numeric>
#>   [1]    chr19 59102001-59102300      * |     -1.00
#>   [2]    chr19 59102301-59102600      * |     -0.75
#>   [3]    chr19 59102601-59102900      * |     -0.50
#>   [4]    chr19 59102901-59103200      * |     -0.25
#>   [5]    chr19 59103201-59103500      * |      0.00
#>   [6]    chr19 59103501-59103800      * |      0.25
#>   [7]    chr17 59103801-59104100      * |      0.50
#>   [8]    chr18 59104101-59104400      * |      0.75
#>   [9]    chr18 59104401-59104700      * |      1.00
#>   -------
#>   seqinfo: 3 sequences from an unspecified genome; no seqlengths
# narrowpeaks
np_file <- system.file("extdata", "demo.narrowPeak.gz",  package="rtracklayer")
gr <- read_narrowpeaks(np_file, genome_info = "hg19")
gr
#> GRanges object with 6 ranges and 6 metadata columns:
#>       seqnames              ranges strand |        name     score signalValue
#>          <Rle>           <IRanges>  <Rle> | <character> <numeric>   <numeric>
#>   [1]    chr19       893203-893659      * |        <NA>         0      16.705
#>   [2]    chr19   45981393-45982404      * |        <NA>         0       7.853
#>   [3]     chr9 137029517-137029957      * |        <NA>         0      19.135
#>   [4]     chr7 148680558-148681129      * |        <NA>         0      10.418
#>   [5]     chr7 148638202-148638769      * |        <NA>         0       8.161
#>   [6]     chr6   52860024-52860589      * |        <NA>         0       9.332
#>          pValue      qValue      peak
#>       <numeric>   <numeric> <integer>
#>   [1]    85.879 3.31125e-82       240
#>   [2]    83.114 1.44562e-79       587
#>   [3]    78.886 1.95499e-75       208
#>   [4]    70.532 3.68097e-67       263
#>   [5]    60.506 2.60538e-57       341
#>   [6]    54.023 7.13040e-51       349
#>   -------
#>   seqinfo: 298 sequences (2 circular) from hg19 genome
```
