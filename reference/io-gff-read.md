# Read a GFF/GTF/GVT file

This is a lightweight wrapper to the import family of functions defined
in rtracklayer.

## Usage

``` r
read_gff(file, col_names = NULL, genome_info = NULL, overlap_ranges = NULL)

read_gff1(file, col_names = NULL, genome_info = NULL, overlap_ranges = NULL)

read_gff2(file, col_names = NULL, genome_info = NULL, overlap_ranges = NULL)

read_gff3(file, col_names = NULL, genome_info = NULL, overlap_ranges = NULL)
```

## Arguments

- file:

  A path to a file or a connection.

- col_names:

  An optional character vector for parsing specific columns in `file`
  that are part of the GFF specification. These should name either fixed
  fields, like source or type, or, for GFF2 and GFF3, any attribute.

- genome_info:

  An optional character string or a Ranges object that contains
  information about the genome build. For example the UCSC identifier
  "hg19" will add build information to the returned GRanges.

- overlap_ranges:

  An optional Ranges object. Only the intervals in the file that overlap
  the Ranges will be returned.

## Value

A GRanges object

a GRanges object

## See also

`rtracklayer::`[`GFFFile()`](https://rdrr.io/pkg/rtracklayer/man/GFFFile-class.html)

## Examples

``` r
test_path <- system.file("tests", package = "rtracklayer")
# gff3
test_gff3 <- file.path(test_path, "genes.gff3")
gr <- read_gff3(test_gff3)
gr
#> GRanges object with 31 ranges and 10 metadata columns:
#>        seqnames      ranges strand |      source     type     score     phase
#>           <Rle>   <IRanges>  <Rle> |    <factor> <factor> <numeric> <integer>
#>    [1]    chr10 92828-95504      - | rtracklayer     gene         5      <NA>
#>    [2]    chr10 92828-95178      - | rtracklayer     mRNA        NA      <NA>
#>    [3]    chr10 92828-95504      - | rtracklayer     mRNA        NA      <NA>
#>    [4]    chr10 92828-94054      - | rtracklayer     exon        NA      <NA>
#>    [5]    chr10 92997-94054      - | rtracklayer     CDS         NA      <NA>
#>    ...      ...         ...    ... .         ...      ...       ...       ...
#>   [27]    chr12 89675-89827      + | rtracklayer     CDS         NA      <NA>
#>   [28]    chr12 90587-90655      + | rtracklayer     exon        NA      <NA>
#>   [29]    chr12 90587-90655      + | rtracklayer     CDS         NA      <NA>
#>   [30]    chr12 90796-91263      + | rtracklayer     exon        NA      <NA>
#>   [31]    chr12 90796-91263      * | rtracklayer     CDS         NA      <NA>
#>                   ID        Name        geneName           Alias      genome
#>          <character> <character>     <character> <CharacterList> <character>
#>    [1] GeneID:347688       TUBB8 tubulin, beta 8  FLJ40100,TUBB8        hg19
#>    [2]           873       TUBB8            <NA>                        <NA>
#>    [3]           872       TUBB8            <NA>                        <NA>
#>    [4]          <NA>        <NA>            <NA>                        <NA>
#>    [5]          <NA>        <NA>            <NA>                        <NA>
#>    ...           ...         ...             ...             ...         ...
#>   [27]          <NA>        <NA>            <NA>                        <NA>
#>   [28]          <NA>        <NA>            <NA>                        <NA>
#>   [29]          <NA>        <NA>            <NA>                        <NA>
#>   [30]          <NA>        <NA>            <NA>                        <NA>
#>   [31]          <NA>        <NA>            <NA>                        <NA>
#>                 Parent
#>        <CharacterList>
#>    [1]                
#>    [2]   GeneID:347688
#>    [3]   GeneID:347688
#>    [4]         872,873
#>    [5]         872,873
#>    ...             ...
#>   [27]            4644
#>   [28]            4644
#>   [29]            4644
#>   [30]            4644
#>   [31]            4644
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
# alternatively with read_gff
gr <- read_gff(test_gff3, genome_info = "hg19")
gr
#> GRanges object with 31 ranges and 10 metadata columns:
#>        seqnames      ranges strand |      source     type     score     phase
#>           <Rle>   <IRanges>  <Rle> |    <factor> <factor> <numeric> <integer>
#>    [1]    chr10 92828-95504      - | rtracklayer     gene         5      <NA>
#>    [2]    chr10 92828-95178      - | rtracklayer     mRNA        NA      <NA>
#>    [3]    chr10 92828-95504      - | rtracklayer     mRNA        NA      <NA>
#>    [4]    chr10 92828-94054      - | rtracklayer     exon        NA      <NA>
#>    [5]    chr10 92997-94054      - | rtracklayer     CDS         NA      <NA>
#>    ...      ...         ...    ... .         ...      ...       ...       ...
#>   [27]    chr12 89675-89827      + | rtracklayer     CDS         NA      <NA>
#>   [28]    chr12 90587-90655      + | rtracklayer     exon        NA      <NA>
#>   [29]    chr12 90587-90655      + | rtracklayer     CDS         NA      <NA>
#>   [30]    chr12 90796-91263      + | rtracklayer     exon        NA      <NA>
#>   [31]    chr12 90796-91263      * | rtracklayer     CDS         NA      <NA>
#>                   ID        Name        geneName           Alias      genome
#>          <character> <character>     <character> <CharacterList> <character>
#>    [1] GeneID:347688       TUBB8 tubulin, beta 8  FLJ40100,TUBB8        hg19
#>    [2]           873       TUBB8            <NA>                        <NA>
#>    [3]           872       TUBB8            <NA>                        <NA>
#>    [4]          <NA>        <NA>            <NA>                        <NA>
#>    [5]          <NA>        <NA>            <NA>                        <NA>
#>    ...           ...         ...             ...             ...         ...
#>   [27]          <NA>        <NA>            <NA>                        <NA>
#>   [28]          <NA>        <NA>            <NA>                        <NA>
#>   [29]          <NA>        <NA>            <NA>                        <NA>
#>   [30]          <NA>        <NA>            <NA>                        <NA>
#>   [31]          <NA>        <NA>            <NA>                        <NA>
#>                 Parent
#>        <CharacterList>
#>    [1]                
#>    [2]   GeneID:347688
#>    [3]   GeneID:347688
#>    [4]         872,873
#>    [5]         872,873
#>    ...             ...
#>   [27]            4644
#>   [28]            4644
#>   [29]            4644
#>   [30]            4644
#>   [31]            4644
#>   -------
#>   seqinfo: 298 sequences (2 circular) from hg19 genome
```
