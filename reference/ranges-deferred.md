# DeferredGenomiRanges objects

Enables deferred reading of files (currently only BAM files) by caching
results after a plyranges verb is called.

## Slots

- `delegate`:

  a GenomicRanges object to be cached

- `ops`:

  A FileOperator object

## See also

[`read_bam()`](https://tidyomics.github.io/plyranges/reference/io-bam-read.md)
