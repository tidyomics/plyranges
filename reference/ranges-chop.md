# Group a GRanges object by introns or gaps

Group a GRanges object by introns or gaps

## Usage

``` r
chop_by_introns(x)

chop_by_gaps(x)
```

## Arguments

- x:

  a GenomicRanges object with a cigar string column

## Value

a GRanges object

## Details

Creates a grouped Ranges object from a cigar string column, for
`chop_by_introns()` will check for the presence of "N" in the cigar
string and create a new column called `intron` where TRUE indicates the
alignment has a skipped region from the reference. For `chop_by_gaps()`
will check for the presence of "N" or "D" in the cigar string and create
a new column called "gaps" where TRUE indicates the alignment has a
deletion from the reference or has an intron.

## Examples

``` r
if (require(pasillaBamSubset)) {
   bamfile <- untreated1_chr4()
   # define a region of interest
   roi <- data.frame(seqnames = "chr4", start = 5e5, end = 7e5) %>%
            as_granges()
   # results in a grouped ranges object
   rng <- read_bam(bamfile) %>% 
            filter_by_overlaps(roi) %>%
            chop_by_gaps()
   # to find ranges that have gaps use filter with `n()`
   rng %>% filter(n() >= 2)
  
}
#> GRanges object with 4258 ranges and 5 metadata columns:
#> Groups: gaps [2129]
#>          seqnames        ranges strand |        which_label       cigar
#>             <Rle>     <IRanges>  <Rle> |              <Rle> <character>
#>      [1]     chr4 501291-501317      - | chr4:500000-700000   27M84N48M
#>      [2]     chr4 501402-501449      - | chr4:500000-700000   27M84N48M
#>      [3]     chr4 501701-501767      - | chr4:500000-700000    67M67N8M
#>      [4]     chr4 501835-501842      - | chr4:500000-700000    67M67N8M
#>      [5]     chr4 501714-501767      - | chr4:500000-700000   54M67N21M
#>      ...      ...           ...    ... .                ...         ...
#>   [4254]     chr4 700351-700408      + | chr4:500000-700000 17M1291N58M
#>   [4255]     chr4 699050-699059      + | chr4:500000-700000 10M1291N65M
#>   [4256]     chr4 700351-700415      + | chr4:500000-700000 10M1291N65M
#>   [4257]     chr4 699050-699059      + | chr4:500000-700000 10M1291N65M
#>   [4258]     chr4 700351-700415      + | chr4:500000-700000 10M1291N65M
#>             qwidth     njunc  gaps
#>          <integer> <integer> <Rle>
#>      [1]        75         1    22
#>      [2]        75         1    22
#>      [3]        75         1    55
#>      [4]        75         1    55
#>      [5]        75         1    56
#>      ...       ...       ...   ...
#>   [4254]        75         1 26613
#>   [4255]        75         1 26627
#>   [4256]        75         1 26627
#>   [4257]        75         1 26629
#>   [4258]        75         1 26629
#>   -------
#>   seqinfo: 8 sequences from an unspecified genome
```
