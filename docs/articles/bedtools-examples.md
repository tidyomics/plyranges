# Using plyranges to perform common genomic data wrangling tasks
{:.no_toc}

<p class="author-name">Stuart Lee</p>

## Contents
{:.no_toc}

-   [Introduction](#introduction)
-   [Data](#data)
-   [Reading BED files](#reading-bed-files)
-   [Intersections and Overlaps](#intersections-and-overlaps)
    -   [Keeping the original features](#keeping-the-original-features)
    -   [Computing the amount of overlap](#computing-the-amount-of-overlap)
    -   [Counting the number of overlaps](#counting-the-number-of-overlaps)
    -   [Excluding by overlaps](#excluding-by-overlaps)
    -   [Excluding by fraction of overlap](#excluding-by-fraction-of-overlap)
-   [Computing Genomic Coverage](#computing-genomic-coverage)
    -   [Coverage Vector](#coverage-vector)
    -   [Coverage histogram](#coverage-histogram)
-   [Composing pipelines](#composing-pipelines)
    -   [Yet another overlaps example](#yet-another-overlaps-example)
{:toc}

## Introduction

Here we introduce the `plyranges` API for manipulating genomic data. The examples follow the [HelloRanges tutorial](http://bioconductor.org/packages/release/bioc/vignettes/HelloRanges/inst/doc/tutorial.pdf). The aim is to show how common analysis tasks that are usually undertaken using `bedtools` can be performed interactively using the `plyranges` API. Throughout this tutorial we will highlight how the `bedtools`, `HelloRanges` and `plyranges` API differ.

The `plyranges` API attempts to be type-consistent, that is if you apply a `plyranges` function to a *Ranges* object it should return a *Ranges* object. This means that intermittant Bioconductor classes such as `Rle` objects are hidden away from the user, and the user can focus on the table-like structure of ordinary *Ranges* objects. A *Ranges* object, in particular a *GRanges* object resembles a BED file and has columns to represent the chromomsome, start, end, and strand information. It can also contain additional columns called metadata containing information about the rows (i.e an exon name or a measurement).

## Data

The data we use for these examples are exported by the `HelloRangesData` package. For more details see that [package's vignette](http://bioconductor.org/packages/release/data/experiment/vignettes/HelloRangesData/inst/doc/intro.pdf). For most examples in this tutorial we will use the RefSeq exons and CpG island annotations from the hg19 genome build.

``` r
exons_bed <- system.file("extdata", "exons.bed", package="HelloRangesData")
cpg_bed <- system.file("extdata", "cpg.bed", package = "HelloRangesData")
```

## Reading BED files

To load the exons file as a *GRanges* object we use the `read_bed` function. We supply the genome build as an additional argument to this function, allowing us to obtain annotation information about each chromosome.

``` r
suppressPackageStartupMessages(library(plyranges))
cpg <- read_bed(cpg_bed, genome_info = "hg19")
cpg
```

    ## GRanges object with 28691 ranges and 1 metadata column:
    ##           seqnames               ranges strand |        name
    ##              <Rle>            <IRanges>  <Rle> | <character>
    ##       [1]     chr1     [ 28736,  29810]      * |    CpG:_116
    ##       [2]     chr1     [135125, 135563]      * |     CpG:_30
    ##       [3]     chr1     [327791, 328229]      * |     CpG:_29
    ##       [4]     chr1     [437152, 438164]      * |     CpG:_84
    ##       [5]     chr1     [449274, 450544]      * |     CpG:_99
    ##       ...      ...                  ...    ... .         ...
    ##   [28687]     chrY [27610116, 27611088]      * |     CpG:_76
    ##   [28688]     chrY [28555536, 28555932]      * |     CpG:_32
    ##   [28689]     chrY [28773316, 28773544]      * |     CpG:_25
    ##   [28690]     chrY [59213795, 59214183]      * |     CpG:_36
    ##   [28691]     chrY [59349267, 59349574]      * |     CpG:_29
    ##   -------
    ##   seqinfo: 93 sequences from hg19 genome

``` r
exons <- read_bed(exons_bed, genome_info = get_genome_info(cpg))
exons
```

    ## GRanges object with 459752 ranges and 2 metadata columns:
    ##            seqnames               ranges strand |
    ##               <Rle>            <IRanges>  <Rle> |
    ##        [1]     chr1       [11874, 12227]      + |
    ##        [2]     chr1       [12613, 12721]      + |
    ##        [3]     chr1       [13221, 14409]      + |
    ##        [4]     chr1       [14362, 14829]      - |
    ##        [5]     chr1       [14970, 15038]      - |
    ##        ...      ...                  ...    ... .
    ##   [459748]     chrY [59338754, 59338859]      + |
    ##   [459749]     chrY [59338754, 59338859]      + |
    ##   [459750]     chrY [59340194, 59340278]      + |
    ##   [459751]     chrY [59342487, 59343488]      + |
    ##   [459752]     chrY [59342487, 59343488]      + |
    ##                                          name     score
    ##                                   <character> <numeric>
    ##        [1]    NR_046018_exon_0_0_chr1_11874_f         0
    ##        [2]    NR_046018_exon_1_0_chr1_12613_f         0
    ##        [3]    NR_046018_exon_2_0_chr1_13221_f         0
    ##        [4]    NR_024540_exon_0_0_chr1_14362_r         0
    ##        [5]    NR_024540_exon_1_0_chr1_14970_r         0
    ##        ...                                ...       ...
    ##   [459748] NM_002186_exon_6_0_chrY_59338754_f         0
    ##   [459749] NM_176786_exon_7_0_chrY_59338754_f         0
    ##   [459750] NM_002186_exon_7_0_chrY_59340194_f         0
    ##   [459751] NM_002186_exon_8_0_chrY_59342487_f         0
    ##   [459752] NM_176786_exon_8_0_chrY_59342487_f         0
    ##   -------
    ##   seqinfo: 93 sequences from hg19 genome

## Intersections and Overlaps

A useful operation on two *Ranges* is to identify where they intersect. By default `bedtools` will return the region of intersection between two genomic tracks. In `plyranges` the intersect operation is a type of overlap join, that is we are finding the common genomic intervals that overlap between the two tracks.

In `plyranges` this is termed an left overlap join:

``` r
intersect_rng <- join_overlap_inner(cpg, exons)
intersect_rng
```

    ## GRanges object with 45500 ranges and 5 metadata columns:
    ##           seqnames               ranges strand |   width.x   width.y
    ##              <Rle>            <IRanges>  <Rle> | <integer> <integer>
    ##       [1]     chr1     [ 29321,  29370]      * |      1075        50
    ##       [2]     chr1     [135125, 135563]      * |       439      4924
    ##       [3]     chr1     [327791, 328229]      * |       439      4143
    ##       [4]     chr1     [327791, 328229]      * |       439      4143
    ##       [5]     chr1     [327791, 328229]      * |       439      1546
    ##       ...      ...                  ...    ... .       ...       ...
    ##   [45496]     chrY [59213949, 59214117]      * |       389       169
    ##   [45497]     chrY [59213949, 59214117]      * |       389       169
    ##   [45498]     chrY [59213949, 59214117]      * |       389       169
    ##   [45499]     chrY [59213949, 59214117]      * |       389       169
    ##   [45500]     chrY [59213949, 59214117]      * |       389       169
    ##                name.x                                name.y     score
    ##           <character>                           <character> <numeric>
    ##       [1]    CpG:_116      NR_024540_exon_10_0_chr1_29321_r         0
    ##       [2]     CpG:_30      NR_039983_exon_0_0_chr1_134773_r         0
    ##       [3]     CpG:_29      NR_028322_exon_2_0_chr1_324439_f         0
    ##       [4]     CpG:_29      NR_028325_exon_2_0_chr1_324439_f         0
    ##       [5]     CpG:_29      NR_028327_exon_3_0_chr1_327036_f         0
    ##       ...         ...                                   ...       ...
    ##   [45496]     CpG:_36 NM_001145149_exon_0_0_chrY_59213949_f         0
    ##   [45497]     CpG:_36 NM_001185183_exon_0_0_chrY_59213949_f         0
    ##   [45498]     CpG:_36    NM_005638_exon_0_0_chrY_59213949_f         0
    ##   [45499]     CpG:_36    NR_033714_exon_0_0_chrY_59213949_f         0
    ##   [45500]     CpG:_36    NR_033715_exon_0_0_chrY_59213949_f         0
    ##   -------
    ##   seqinfo: 93 sequences from hg19 genome

In the `HelloRanges` API (we also see the `bedtools` code here) the equivalent operation is:

``` r
suppressPackageStartupMessages(library(HelloRanges))
code <- bedtools_intersect("-a cpg.bed -b exons.bed -g hg19")
code
```

    ## {
    ##     genome <- Seqinfo(genome = "hg19")
    ##     gr_a <- import("cpg.bed", genome = genome)
    ##     gr_b <- import("exons.bed", genome = genome)
    ##     pairs <- findOverlapPairs(gr_a, gr_b, ignore.strand = TRUE)
    ##     ans <- pintersect(pairs, ignore.strand = TRUE)
    ##     ans
    ## }

### Keeping the original features

By default an overlap join will keep the information in both the query and subject ranges.

### Computing the amount of overlap

To compute the amount of overlap we store the width of the intersecting ranges as an additional column using the `mutate` operator.

``` r
intersect_rng <- intersect_rng %>% 
  mutate(overlap_width = width)
```

The equivalent bedtools/HelloRanges command is:

``` r
bedtools_intersect("-a cpg.bed -b exons.bed -g hg19 -wo")
```

    ## {
    ##     genome <- Seqinfo(genome = "hg19")
    ##     gr_a <- import("cpg.bed", genome = genome)
    ##     gr_b <- import("exons.bed", genome = genome)
    ##     pairs <- findOverlapPairs(gr_a, gr_b, ignore.strand = TRUE)
    ##     ans <- pairs
    ##     mcols(ans)$overlap_width <- width(pintersect(ans, ignore.strand = TRUE))
    ##     ans
    ## }

### Counting the number of overlaps

To add the count of the number of times each cpg island overlaps an exon we can use `mutate` with the `count_overlaps` function.

``` r
cpg %>% 
  mutate(overlaps = count_overlaps(., exons))
```

    ## GRanges object with 28691 ranges and 2 metadata columns:
    ##           seqnames               ranges strand |        name  overlaps
    ##              <Rle>            <IRanges>  <Rle> | <character> <integer>
    ##       [1]     chr1     [ 28736,  29810]      * |    CpG:_116         1
    ##       [2]     chr1     [135125, 135563]      * |     CpG:_30         1
    ##       [3]     chr1     [327791, 328229]      * |     CpG:_29         3
    ##       [4]     chr1     [437152, 438164]      * |     CpG:_84         0
    ##       [5]     chr1     [449274, 450544]      * |     CpG:_99         0
    ##       ...      ...                  ...    ... .         ...       ...
    ##   [28687]     chrY [27610116, 27611088]      * |     CpG:_76         0
    ##   [28688]     chrY [28555536, 28555932]      * |     CpG:_32         0
    ##   [28689]     chrY [28773316, 28773544]      * |     CpG:_25         0
    ##   [28690]     chrY [59213795, 59214183]      * |     CpG:_36         5
    ##   [28691]     chrY [59349267, 59349574]      * |     CpG:_29         0
    ##   -------
    ##   seqinfo: 93 sequences from hg19 genome

The equivalent bedtools operation is

``` r
bedtools_intersect("-a cpg.bed -b exons.bed -g hg19 -c")
```

    ## {
    ##     genome <- Seqinfo(genome = "hg19")
    ##     gr_a <- import("cpg.bed", genome = genome)
    ##     gr_b <- import("exons.bed", genome = genome)
    ##     ans <- gr_a
    ##     mcols(ans)$overlap_count <- countOverlaps(gr_a, gr_b, ignore.strand = TRUE)
    ##     ans
    ## }

### Excluding by overlaps

We can use `filter` along with the infix operator `%over` to only include CpG islands that do not overlap any exons.

``` r
cpg %>% 
  filter_by_non_overlaps(., exons)
```

    ## GRanges object with 9846 ranges and 1 metadata column:
    ##          seqnames               ranges strand |        name
    ##             <Rle>            <IRanges>  <Rle> | <character>
    ##      [1]     chr1     [437152, 438164]      * |     CpG:_84
    ##      [2]     chr1     [449274, 450544]      * |     CpG:_99
    ##      [3]     chr1     [533220, 534114]      * |     CpG:_94
    ##      [4]     chr1     [544739, 546649]      * |    CpG:_171
    ##      [5]     chr1     [801976, 802338]      * |     CpG:_24
    ##      ...      ...                  ...    ... .         ...
    ##   [9842]     chrY [26351344, 26352316]      * |     CpG:_76
    ##   [9843]     chrY [27610116, 27611088]      * |     CpG:_76
    ##   [9844]     chrY [28555536, 28555932]      * |     CpG:_32
    ##   [9845]     chrY [28773316, 28773544]      * |     CpG:_25
    ##   [9846]     chrY [59349267, 59349574]      * |     CpG:_29
    ##   -------
    ##   seqinfo: 93 sequences from hg19 genome

The equivalent bedtools operation is

``` r
bedtools_intersect("-a cpg.bed -b exons.bed -g hg19 -v")
```

    ## {
    ##     genome <- Seqinfo(genome = "hg19")
    ##     gr_a <- import("cpg.bed", genome = genome)
    ##     gr_b <- import("exons.bed", genome = genome)
    ##     subsetByOverlaps(gr_a, gr_b, invert = TRUE, ignore.strand = TRUE)
    ## }

### Excluding by fraction of overlap

A common operation is to filter ranges by the fraction they overlap a given query or subject. We can achieve this with an overlap inner join and mutate. Since an overlap inner join will keep the width of the original query ranges we can combine it with filter. In this example we filter ranges where the proportion of overlap on cpg islands is less than half.

``` r
olap <- join_overlap_inner(cpg, exons) %>% 
  filter(width / width.x  >= 0.5)
olap
```

    ## GRanges object with 12669 ranges and 5 metadata columns:
    ##           seqnames               ranges strand |   width.x   width.y
    ##              <Rle>            <IRanges>  <Rle> | <integer> <integer>
    ##       [1]     chr1     [135125, 135563]      * |       439      4924
    ##       [2]     chr1     [327791, 328229]      * |       439      4143
    ##       [3]     chr1     [327791, 328229]      * |       439      4143
    ##       [4]     chr1     [327791, 328229]      * |       439      1546
    ##       [5]     chr1     [788864, 789211]      * |       348      6056
    ##       ...      ...                  ...    ... .       ...       ...
    ##   [12665]     chrY [26979967, 26980116]      * |       227       310
    ##   [12666]     chrY [26979967, 26980116]      * |       227       310
    ##   [12667]     chrY [26979967, 26980116]      * |       227       310
    ##   [12668]     chrY [26979990, 26980116]      * |       227       287
    ##   [12669]     chrY [26979990, 26980116]      * |       227       287
    ##                name.x                                name.y     score
    ##           <character>                           <character> <numeric>
    ##       [1]     CpG:_30      NR_039983_exon_0_0_chr1_134773_r         0
    ##       [2]     CpG:_29      NR_028322_exon_2_0_chr1_324439_f         0
    ##       [3]     CpG:_29      NR_028325_exon_2_0_chr1_324439_f         0
    ##       [4]     CpG:_29      NR_028327_exon_3_0_chr1_327036_f         0
    ##       [5]     CpG:_28      NR_047519_exon_5_0_chr1_788771_f         0
    ##       ...         ...                                   ...       ...
    ##   [12665]     CpG:_21 NM_001005375_exon_0_0_chrY_26979967_f         0
    ##   [12666]     CpG:_21    NM_020364_exon_0_0_chrY_26979967_f         0
    ##   [12667]     CpG:_21    NM_020420_exon_0_0_chrY_26979967_f         0
    ##   [12668]     CpG:_21 NM_001005785_exon_0_0_chrY_26979990_f         0
    ##   [12669]     CpG:_21 NM_001005786_exon_0_0_chrY_26979990_f         0
    ##   -------
    ##   seqinfo: 93 sequences from hg19 genome

The equivalent bedtools code is:

``` r
bedtools_intersect("-a cpg.bed -b exons.bed -g hg19 -f 0.5 -wo")
```

    ## {
    ##     genome <- Seqinfo(genome = "hg19")
    ##     gr_a <- import("cpg.bed", genome = genome)
    ##     gr_b <- import("exons.bed", genome = genome)
    ##     pairs <- findOverlapPairs(gr_a, gr_b, ignore.strand = TRUE)
    ##     olap <- pintersect(pairs, ignore.strand = TRUE)
    ##     keep <- width(olap)/width(first(pairs)) >= 0.5
    ##     pairs <- pairs[keep]
    ##     ans <- pairs
    ##     mcols(ans)$overlap_width <- width(olap)[keep]
    ##     ans
    ## }

## Computing Genomic Coverage

Often we are interested in counting the number of features over an entire genome that overlap each other. In `plyranges` this is done with the `set_coverage` function in combination with other functions to manipulate the results. This function will always return a new *Ranges* object with a column called score. Any other variables associated with the input *Ranges* will be dropped.

### Coverage Vector

By default, `set_coverage` is equivalent to `bedtools genomecov -bga` and hence will return scores equal to zero. Filtering all scores above zero is equivalent to `bedtools genomecov -bg`:

``` r
cvg <- exons %>% 
  set_coverage() %>% 
  filter(score > 0)
cvg
```

    ## GRanges object with 245860 ranges and 1 metadata column:
    ##                  seqnames           ranges strand |     score
    ##                     <Rle>        <IRanges>  <Rle> | <integer>
    ##        [1]           chr1   [11874, 12227]      * |         1
    ##        [2]           chr1   [12613, 12721]      * |         1
    ##        [3]           chr1   [13221, 14361]      * |         1
    ##        [4]           chr1   [14362, 14409]      * |         2
    ##        [5]           chr1   [14410, 14829]      * |         1
    ##        ...            ...              ...    ... .       ...
    ##   [245856] chrUn_gl000228 [109202, 110581]      * |         6
    ##   [245857] chrUn_gl000228 [112508, 112604]      * |         6
    ##   [245858] chrUn_gl000228 [112605, 113887]      * |         7
    ##   [245859] chrUn_gl000228 [114024, 114115]      * |         1
    ##   [245860] chrUn_gl000228 [114478, 114676]      * |         1
    ##   -------
    ##   seqinfo: 93 sequences from an unspecified genome

``` r
bedtools_genomecov("-i exons.bed -g hg19.genome -bg")
```

    ## {
    ##     genome <- import("hg19.genome")
    ##     gr_a <- import("exons.bed", genome = genome)
    ##     cov <- coverage(gr_a)
    ##     ans <- GRanges(cov)
    ##     ans <- subset(ans, score > 0)
    ##     ans
    ## }

### Coverage histogram

We can construct a histogram over all seqeunces in a genome using the `group_by` and `summarise` operations along with `left_join`.

First, we compute the coverage over the exons and construct a `tibble` from the genome build information using `select` with `.drop_ranges = FALSE`.

``` r
suppressPackageStartupMessages(library(dplyr))
cvg <- exons %>% set_coverage()

## convert the sequence annotation to a tibble
hg19 <- get_genome_info(cvg) %>%  
 select(seqnames, width, .drop_ranges = TRUE)
```

Then to count the number of bases corresponding to a score we sum over the width of each range in each chromosome. Then we join the resulting `tibble` to the annotation `tibble` called `hg19`. Note that in the sums we have coerced integer variables to doubles to avoid overflow.

``` r
cvg_hist_by_seq <- cvg %>%
  group_by(seqnames, score) %>%
  summarise(count = sum(as.numeric(width))) %>%  
  left_join(., hg19, by = "seqnames") %>% 
  mutate(fraction = count / width)

cvg_hist_by_seq
```

    ## # A tibble: 608 x 5
    ##    seqnames score     count     width  fraction
    ##      <fctr> <int>     <dbl>     <int>     <dbl>
    ##  1     chr1     0 241996316 249250621 0.9708955
    ##  2     chr2     0 238060317 243199373 0.9788690
    ##  3     chr3     0 193800642 198022430 0.9786803
    ##  4     chr4     0 188206114 191154276 0.9845771
    ##  5     chr5     0 177466741 180915260 0.9809385
    ##  6     chr6     0 167349524 171115067 0.9779941
    ##  7     chr7     0 155623839 159138663 0.9779135
    ##  8     chr8     0 143743330 146364022 0.9820947
    ##  9     chr9     0 138299162 141213431 0.9793627
    ## 10    chr10     0 132497709 135534747 0.9775922
    ## # ... with 598 more rows

Similarly for the genome-wide coverage histogram, we perform the same operation but do not group over `seqnames`. Finally we bind the resulting `tibbles` together.

``` r
cvg_hist_total <- cvg %>%
  group_by(score) %>% 
  summarise(count = sum(as.numeric(width))) %>% 
  mutate(seqnames = "genome") %>%
  left_join(.,
            hg19 %>% 
              summarise(width = sum(as.numeric(width))) %>% 
              mutate(seqnames = "genome")) %>%
  mutate(fraction = count / width)
```

    ## Joining, by = "seqnames"

``` r
cvg_hist_total
```

    ## # A tibble: 45 x 5
    ##    score      count seqnames      width     fraction
    ##    <int>      <dbl>    <chr>      <dbl>        <dbl>
    ##  1     0 3062406951   genome 3137161264 0.9761713515
    ##  2     1   44120515   genome 3137161264 0.0140638339
    ##  3     2   15076446   genome 3137161264 0.0048057606
    ##  4     3    7294047   genome 3137161264 0.0023250469
    ##  5     4    3650324   genome 3137161264 0.0011635755
    ##  6     5    1926397   genome 3137161264 0.0006140574
    ##  7     6    1182623   genome 3137161264 0.0003769723
    ##  8     7     574102   genome 3137161264 0.0001830005
    ##  9     8     353352   genome 3137161264 0.0001126343
    ## 10     9     152653   genome 3137161264 0.0000486596
    ## # ... with 35 more rows

``` r
cvg_hist <- bind_rows(cvg_hist_by_seq, cvg_hist_total)
cvg_hist
```

    ## # A tibble: 653 x 5
    ##    seqnames score     count     width  fraction
    ##       <chr> <int>     <dbl>     <dbl>     <dbl>
    ##  1     chr1     0 241996316 249250621 0.9708955
    ##  2     chr2     0 238060317 243199373 0.9788690
    ##  3     chr3     0 193800642 198022430 0.9786803
    ##  4     chr4     0 188206114 191154276 0.9845771
    ##  5     chr5     0 177466741 180915260 0.9809385
    ##  6     chr6     0 167349524 171115067 0.9779941
    ##  7     chr7     0 155623839 159138663 0.9779135
    ##  8     chr8     0 143743330 146364022 0.9820947
    ##  9     chr9     0 138299162 141213431 0.9793627
    ## 10    chr10     0 132497709 135534747 0.9775922
    ## # ... with 643 more rows

Although this is slightly more verbose than the `bedtools` or `HelloRanges` approach the `plyranges` code makes the actions being performed on the input `Ranges` explicit:

``` r
bedtools_genomecov("-i exons.bed -g hg19.genome")
```

    ## {
    ##     genome <- import("hg19.genome")
    ##     gr_a <- import("exons.bed", genome = genome)
    ##     cov <- coverage(gr_a)
    ##     tablist <- List(lapply(cov, table))
    ##     mcols(tablist)$len <- lengths(cov, use.names = FALSE)
    ##     covhist <- stack(tablist, "seqnames", "count", "coverage")
    ##     margin <- aggregate(covhist, ~coverage, count = sum(NumericList(count)))[-1L]
    ##     margin <- DataFrame(seqnames = Rle("genome"), margin, len = sum(as.numeric(lengths(cov))))
    ##     covhist <- rbind(covhist, margin)
    ##     ans <- within(covhist, fraction <- count/len)
    ##     ans
    ## }

## Composing pipelines

### Yet another overlaps example

Here we perform another example where we chain a filter and then find overlapping ranges between two ranges. In this example we reduce the exons to have zero and then find the ranges that overlap between the filtered exon ranges and the cpg islands.

``` r
overlaps <- exons %>% 
  filter(score == 0L) %>% 
  join_overlap_inner(.,  cpg) 

overlaps
```

    ## GRanges object with 45500 ranges and 5 metadata columns:
    ##           seqnames               ranges strand |   width.x   width.y
    ##              <Rle>            <IRanges>  <Rle> | <integer> <integer>
    ##       [1]     chr1     [ 29321,  29370]      - |        50      1075
    ##       [2]     chr1     [135125, 135563]      - |      4924       439
    ##       [3]     chr1     [327791, 328229]      + |      4143       439
    ##       [4]     chr1     [327791, 328229]      + |      4143       439
    ##       [5]     chr1     [327791, 328229]      + |      1546       439
    ##       ...      ...                  ...    ... .       ...       ...
    ##   [45496]     chrY [59213949, 59214117]      + |       169       389
    ##   [45497]     chrY [59213949, 59214117]      + |       169       389
    ##   [45498]     chrY [59213949, 59214117]      + |       169       389
    ##   [45499]     chrY [59213949, 59214117]      + |       169       389
    ##   [45500]     chrY [59213949, 59214117]      + |       169       389
    ##                                          name.x     score      name.y
    ##                                     <character> <numeric> <character>
    ##       [1]      NR_024540_exon_10_0_chr1_29321_r         0    CpG:_116
    ##       [2]      NR_039983_exon_0_0_chr1_134773_r         0     CpG:_30
    ##       [3]      NR_028322_exon_2_0_chr1_324439_f         0     CpG:_29
    ##       [4]      NR_028325_exon_2_0_chr1_324439_f         0     CpG:_29
    ##       [5]      NR_028327_exon_3_0_chr1_327036_f         0     CpG:_29
    ##       ...                                   ...       ...         ...
    ##   [45496] NM_001145149_exon_0_0_chrY_59213949_f         0     CpG:_36
    ##   [45497] NM_001185183_exon_0_0_chrY_59213949_f         0     CpG:_36
    ##   [45498]    NM_005638_exon_0_0_chrY_59213949_f         0     CpG:_36
    ##   [45499]    NR_033714_exon_0_0_chrY_59213949_f         0     CpG:_36
    ##   [45500]    NR_033715_exon_0_0_chrY_59213949_f         0     CpG:_36
    ##   -------
    ##   seqinfo: 93 sequences from hg19 genome

We can also compute the coverage histogram of exons over cpg islands and then plot results as an ecdf.

``` r
cvg_over_exons <- exons %>% 
  set_coverage() %>%
  join_overlap_inner(., cpg) %>%
  mutate(total = sum(as.numeric(width))) %>% 
  group_by(score) %>%
  summarise(count = sum(as.numeric(width)),
            fraction = count / unique(total)) 

cvg_over_exons
```

    ## # A tibble: 28 x 3
    ##    score    count     fraction
    ##    <int>    <dbl>        <dbl>
    ##  1     0 14607020 0.6687356377
    ##  2     1  4644079 0.2126142862
    ##  3     2  1431318 0.0655283114
    ##  4     3   581001 0.0265992704
    ##  5     4   245183 0.0112249186
    ##  6     5   137255 0.0062837807
    ##  7     6    84053 0.0038480975
    ##  8     7    41162 0.0018844704
    ##  9     8    32327 0.0014799882
    ## 10     9     9186 0.0004205516
    ## # ... with 18 more rows

``` r
library(ggplot2)
ggplot(cvg_over_exons, aes(x = score, y = 1 - cumsum(fraction))) +
  geom_step() +
  xlim(c(0,25))  +
  xlab("coverage") +
  ylab("fraction of bp > coverage")
```

![](bedtools-examples_files/figure-markdown_github/unnamed-chunk-15-1.png)
