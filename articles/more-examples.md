# Additional examples of plyranges

## Quick overview

### About `Ranges`

`Ranges` objects can either represent sets of integers as `IRanges`
(which have start, end and width attributes) or represent genomic
intervals (which have additional attributes, sequence name, and strand)
as `GRanges`. In addition, both types of `Ranges` can store information
about their intervals as metadata columns (for example GC content over a
genomic interval).

`Ranges` objects follow the tidy data principle: each row of a `Ranges`
object corresponds to an interval, while each column will represent a
variable about that interval, and generally each object will represent a
single unit of observation (like gene annotations).

We can construct a `IRanges` object from a `data.frame` with a `start`
or `width` using the
[`as_iranges()`](https://tidyomics.github.io/plyranges/reference/ranges-construct.md)
method.

``` r

library(plyranges)
df <- data.frame(start = 1:5, width = 5)
as_iranges(df)
```

    ## IRanges object with 5 ranges and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]         1         5         5
    ##   [2]         2         6         5
    ##   [3]         3         7         5
    ##   [4]         4         8         5
    ##   [5]         5         9         5

``` r

# alternatively with end
df <- data.frame(start = 1:5, end = 5:9)
as_iranges(df)
```

    ## IRanges object with 5 ranges and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]         1         5         5
    ##   [2]         2         6         5
    ##   [3]         3         7         5
    ##   [4]         4         8         5
    ##   [5]         5         9         5

We can also construct a `GRanges` object in a similar manner. Note that
a `GRanges` object requires at least a seqnames column to be present in
the data.frame (but not necessarily a strand column).

``` r

df <- data.frame(seqnames = c("chr1", "chr2", "chr2", "chr1", "chr2"),
                 start = 1:5,
                 width = 5)
as_granges(df)
```

    ## GRanges object with 5 ranges and 0 metadata columns:
    ##       seqnames    ranges strand
    ##          <Rle> <IRanges>  <Rle>
    ##   [1]     chr1       1-5      *
    ##   [2]     chr2       2-6      *
    ##   [3]     chr2       3-7      *
    ##   [4]     chr1       4-8      *
    ##   [5]     chr2       5-9      *
    ##   -------
    ##   seqinfo: 2 sequences from an unspecified genome; no seqlengths

``` r

# strand can be specified with `+`, `*` (mising) and `-`
df$strand <- c("+", "+", "-", "-", "*")
as_granges(df)
```

    ## GRanges object with 5 ranges and 0 metadata columns:
    ##       seqnames    ranges strand
    ##          <Rle> <IRanges>  <Rle>
    ##   [1]     chr1       1-5      +
    ##   [2]     chr2       2-6      +
    ##   [3]     chr2       3-7      -
    ##   [4]     chr1       4-8      -
    ##   [5]     chr2       5-9      *
    ##   -------
    ##   seqinfo: 2 sequences from an unspecified genome; no seqlengths

## Example: finding GWAS hits that overlap known exons

Let’s look at a more a realistic example (taken from HelloRanges
vignette).

Suppose we have two *GRanges* objects: one containing coordinates of
known exons and another containing SNPs from a GWAS.

The first and last 5 exons are printed below, there are two additional
columns corresponding to the exon name, and a score.

We could check the number of exons per chromosome using `group_by` and
`summarise`.

``` r

exons
```

    ## GRanges object with 459752 ranges and 2 metadata columns:
    ##            seqnames            ranges strand |                   name     score
    ##               <Rle>         <IRanges>  <Rle> |            <character> <numeric>
    ##        [1]     chr1       11874-12227      + | NR_046018_exon_0_0_c..         0
    ##        [2]     chr1       12613-12721      + | NR_046018_exon_1_0_c..         0
    ##        [3]     chr1       13221-14409      + | NR_046018_exon_2_0_c..         0
    ##        [4]     chr1       14362-14829      - | NR_024540_exon_0_0_c..         0
    ##        [5]     chr1       14970-15038      - | NR_024540_exon_1_0_c..         0
    ##        ...      ...               ...    ... .                    ...       ...
    ##   [459748]     chrY 59338754-59338859      + | NM_002186_exon_6_0_c..         0
    ##   [459749]     chrY 59338754-59338859      + | NM_176786_exon_7_0_c..         0
    ##   [459750]     chrY 59340194-59340278      + | NM_002186_exon_7_0_c..         0
    ##   [459751]     chrY 59342487-59343488      + | NM_002186_exon_8_0_c..         0
    ##   [459752]     chrY 59342487-59343488      + | NM_176786_exon_8_0_c..         0
    ##   -------
    ##   seqinfo: 93 sequences from an unspecified genome; no seqlengths

``` r

exons %>%
  group_by(seqnames) %>%
  summarise(n = n())
```

    ## DataFrame with 49 rows and 2 columns
    ##           seqnames         n
    ##              <Rle> <integer>
    ## 1             chr1     43366
    ## 2            chr10     19420
    ## 3            chr11     24476
    ## 4            chr12     24949
    ## 5            chr13      7974
    ## ...            ...       ...
    ## 45  chrUn_gl000222        20
    ## 46  chrUn_gl000223        22
    ## 47  chrUn_gl000228        85
    ## 48            chrX     18173
    ## 49            chrY      4128

Next we create a column representing the transcript_id with `mutate`:

``` r

exons <- exons %>%
  mutate(tx_id = sub("_exon.*", "", name))
```

To find all GWAS SNPs that overlap exons, we use `join_overlap_inner`.
This will create a new *GRanges* with the coordinates of SNPs that
overlap exons, as well as metadata from both objects.

``` r

olap <- join_overlap_inner(gwas, exons)
olap
```

    ## GRanges object with 3439 ranges and 4 metadata columns:
    ##          seqnames    ranges strand |      name.x                 name.y
    ##             <Rle> <IRanges>  <Rle> | <character>            <character>
    ##      [1]     chr1   1079198      * |  rs11260603 NR_038869_exon_2_0_c..
    ##      [2]     chr1   1247494      * |     rs12103 NM_001256456_exon_1_..
    ##      [3]     chr1   1247494      * |     rs12103 NM_001256460_exon_1_..
    ##      [4]     chr1   1247494      * |     rs12103 NM_001256462_exon_1_..
    ##      [5]     chr1   1247494      * |     rs12103 NM_001256463_exon_1_..
    ##      ...      ...       ...    ... .         ...                    ...
    ##   [3435]     chrX 153764217      * |   rs1050828 NM_001042351_exon_9_..
    ##   [3436]     chrX 153764217      * |   rs1050828 NM_000402_exon_9_0_c..
    ##   [3437]     chrX 153764217      * |   rs1050828 NM_001042351_exon_9_..
    ##   [3438]     chrX 153764217      * |   rs1050828 NM_000402_exon_9_0_c..
    ##   [3439]     chrX 153764217      * |   rs1050828 NM_001042351_exon_9_..
    ##              score        tx_id
    ##          <numeric>  <character>
    ##      [1]         0    NR_038869
    ##      [2]         0 NM_001256456
    ##      [3]         0 NM_001256460
    ##      [4]         0 NM_001256462
    ##      [5]         0 NM_001256463
    ##      ...       ...          ...
    ##   [3435]         0 NM_001042351
    ##   [3436]         0    NM_000402
    ##   [3437]         0 NM_001042351
    ##   [3438]         0    NM_000402
    ##   [3439]         0 NM_001042351
    ##   -------
    ##   seqinfo: 93 sequences from an unspecified genome; no seqlengths

For each SNP we can count the number of times it overlaps a transcript.

``` r

olap %>%
  group_by(name.x, tx_id) %>%
  summarise(n = n())
```

    ## DataFrame with 1619 rows and 3 columns
    ##           name.x        tx_id         n
    ##      <character>  <character> <integer>
    ## 1     rs10043775 NM_001271723         1
    ## 2     rs10043775    NM_030793         1
    ## 3        rs10078 NM_001242412         1
    ## 4        rs10078    NM_020731         1
    ## 5        rs10089    NM_001046         1
    ## ...          ...          ...       ...
    ## 1615   rs9906595 NM_001008777         1
    ## 1616      rs9948    NM_017623         1
    ## 1617      rs9948    NM_199078         1
    ## 1618    rs995030    NM_000899         4
    ## 1619    rs995030    NM_003994         4

We can also generate 2bp splice sites on either side of the exon using
`flank_left` and `flank_right`. We add a column indicating the side of
flanking for illustrative purposes. The `interweave` function pairs the
left and right ranges objects.

``` r

left_ss <- flank_left(exons, 2L)
right_ss <- flank_right(exons, 2L)
all_ss <- interweave(left_ss, right_ss, .id = "side")
all_ss
```

    ## GRanges object with 919504 ranges and 4 metadata columns:
    ##            seqnames            ranges strand |                   name     score
    ##               <Rle>         <IRanges>  <Rle> |            <character> <numeric>
    ##        [1]     chr1       11872-11873      + | NR_046018_exon_0_0_c..         0
    ##        [2]     chr1       12228-12229      + | NR_046018_exon_0_0_c..         0
    ##        [3]     chr1       12611-12612      + | NR_046018_exon_1_0_c..         0
    ##        [4]     chr1       12722-12723      + | NR_046018_exon_1_0_c..         0
    ##        [5]     chr1       13219-13220      + | NR_046018_exon_2_0_c..         0
    ##        ...      ...               ...    ... .                    ...       ...
    ##   [919500]     chrY 59340279-59340280      + | NM_002186_exon_7_0_c..         0
    ##   [919501]     chrY 59342485-59342486      + | NM_002186_exon_8_0_c..         0
    ##   [919502]     chrY 59343489-59343490      + | NM_002186_exon_8_0_c..         0
    ##   [919503]     chrY 59342485-59342486      + | NM_176786_exon_8_0_c..         0
    ##   [919504]     chrY 59343489-59343490      + | NM_176786_exon_8_0_c..         0
    ##                  tx_id        side
    ##            <character> <character>
    ##        [1]   NR_046018        left
    ##        [2]   NR_046018       right
    ##        [3]   NR_046018        left
    ##        [4]   NR_046018       right
    ##        [5]   NR_046018        left
    ##        ...         ...         ...
    ##   [919500]   NM_002186       right
    ##   [919501]   NM_002186        left
    ##   [919502]   NM_002186       right
    ##   [919503]   NM_176786        left
    ##   [919504]   NM_176786       right
    ##   -------
    ##   seqinfo: 93 sequences from an unspecified genome; no seqlengths

## Session information

``` r

sessionInfo()
```

    ## R version 4.5.2 (2025-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] plyranges_1.31.7     dplyr_1.2.1          GenomicRanges_1.62.1
    ## [4] Seqinfo_1.0.0        IRanges_2.44.0       S4Vectors_0.48.1    
    ## [7] BiocGenerics_0.56.0  generics_0.1.4       BiocStyle_2.38.0    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] SummarizedExperiment_1.40.0 rjson_0.2.23               
    ##  [3] xfun_0.58                   bslib_0.11.0               
    ##  [5] htmlwidgets_1.6.4           Biobase_2.70.0             
    ##  [7] lattice_0.22-9              vctrs_0.7.3                
    ##  [9] tools_4.5.2                 bitops_1.0-9               
    ## [11] curl_7.1.0                  parallel_4.5.2             
    ## [13] tibble_3.3.1                pkgconfig_2.0.3            
    ## [15] Matrix_1.7-5                desc_1.4.3                 
    ## [17] cigarillo_1.0.0             lifecycle_1.0.5            
    ## [19] compiler_4.5.2              Rsamtools_2.26.0           
    ## [21] textshaping_1.0.5           Biostrings_2.78.0          
    ## [23] codetools_0.2-20            htmltools_0.5.9            
    ## [25] sass_0.4.10                 RCurl_1.98-1.18            
    ## [27] yaml_2.3.12                 pillar_1.11.1              
    ## [29] pkgdown_2.2.0               crayon_1.5.3               
    ## [31] jquerylib_0.1.4             BiocParallel_1.44.0        
    ## [33] DelayedArray_0.36.1         cachem_1.1.0               
    ## [35] abind_1.4-8                 tidyselect_1.2.1           
    ## [37] digest_0.6.39               restfulr_0.0.16            
    ## [39] bookdown_0.46               fastmap_1.2.0              
    ## [41] grid_4.5.2                  cli_3.6.6                  
    ## [43] SparseArray_1.10.10         magrittr_2.0.5             
    ## [45] S4Arrays_1.10.1             XML_3.99-0.23              
    ## [47] withr_3.0.2                 rmarkdown_2.31             
    ## [49] XVector_0.50.0              httr_1.4.8                 
    ## [51] matrixStats_1.5.0           otel_0.2.0                 
    ## [53] ragg_1.5.2                  evaluate_1.0.5             
    ## [55] knitr_1.51                  BiocIO_1.20.0              
    ## [57] rtracklayer_1.70.1          rlang_1.2.0                
    ## [59] glue_1.8.1                  BiocManager_1.30.27        
    ## [61] jsonlite_2.0.0              R6_2.6.1                   
    ## [63] MatrixGenerics_1.22.0       GenomicAlignments_1.46.0   
    ## [65] systemfonts_1.3.2           fs_2.1.0
