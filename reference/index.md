# Package index

## About plyranges

- [`plyranges`](https://tidyomics.github.io/plyranges/reference/plyranges-package.md)
  [`plyranges-package`](https://tidyomics.github.io/plyranges/reference/plyranges-package.md)
  : plyranges: a grammar of genomic data manipulation

## Ranges

- [`as_iranges()`](https://tidyomics.github.io/plyranges/reference/ranges-construct.md)
  [`as_granges()`](https://tidyomics.github.io/plyranges/reference/ranges-construct.md)
  : Construct a I/GRanges object from a tibble or data.frame
- [`interweave()`](https://tidyomics.github.io/plyranges/reference/ranges-interweave.md)
  : Interweave a pair of Ranges objects together
- [`bind_ranges()`](https://tidyomics.github.io/plyranges/reference/ranges-bind.md)
  : Combine Ranges by concatentating them together
- [`tile_ranges()`](https://tidyomics.github.io/plyranges/reference/ranges-tile.md)
  [`slide_ranges()`](https://tidyomics.github.io/plyranges/reference/ranges-tile.md)
  : Slide or tile over a Ranges object

## Arithmetic

- [`anchor()`](https://tidyomics.github.io/plyranges/reference/ranges-anchor.md)
  [`unanchor()`](https://tidyomics.github.io/plyranges/reference/ranges-anchor.md)
  [`anchor_start()`](https://tidyomics.github.io/plyranges/reference/ranges-anchor.md)
  [`anchor_end()`](https://tidyomics.github.io/plyranges/reference/ranges-anchor.md)
  [`anchor_center()`](https://tidyomics.github.io/plyranges/reference/ranges-anchor.md)
  [`anchor_centre()`](https://tidyomics.github.io/plyranges/reference/ranges-anchor.md)
  [`anchor_3p()`](https://tidyomics.github.io/plyranges/reference/ranges-anchor.md)
  [`anchor_5p()`](https://tidyomics.github.io/plyranges/reference/ranges-anchor.md)
  : Anchored Ranges objects
- [`stretch()`](https://tidyomics.github.io/plyranges/reference/stretch.md)
  : Stretch a genomic interval
- [`shift_left()`](https://tidyomics.github.io/plyranges/reference/ranges-shift.md)
  [`shift_right()`](https://tidyomics.github.io/plyranges/reference/ranges-shift.md)
  [`shift_upstream()`](https://tidyomics.github.io/plyranges/reference/ranges-shift.md)
  [`shift_downstream()`](https://tidyomics.github.io/plyranges/reference/ranges-shift.md)
  : Shift all coordinates in a genomic interval left or right, upstream
  or downstream
- [`flank_left()`](https://tidyomics.github.io/plyranges/reference/ranges-flank.md)
  [`flank_right()`](https://tidyomics.github.io/plyranges/reference/ranges-flank.md)
  [`flank_upstream()`](https://tidyomics.github.io/plyranges/reference/ranges-flank.md)
  [`flank_downstream()`](https://tidyomics.github.io/plyranges/reference/ranges-flank.md)
  : Generate flanking regions
- [`chop_by_introns()`](https://tidyomics.github.io/plyranges/reference/ranges-chop.md)
  [`chop_by_gaps()`](https://tidyomics.github.io/plyranges/reference/ranges-chop.md)
  : Group a GRanges object by introns or gaps
- [`n()`](https://tidyomics.github.io/plyranges/reference/n.md) :
  Compute the number of ranges in each group.
- [`n_distinct()`](https://tidyomics.github.io/plyranges/reference/n_distinct.md)
  : Compute the number of distinct unique values in a vector or List

## Core verbs

- [`mutate(`*`<Ranges>`*`)`](https://tidyomics.github.io/plyranges/reference/mutate-ranges.md)
  : Modify a Ranges object

- [`summarise(`*`<Ranges>`*`)`](https://tidyomics.github.io/plyranges/reference/ranges-summarise.md)
  : Reduce multiple values in a Ranges down to a single value

- [`filter(`*`<Ranges>`*`)`](https://tidyomics.github.io/plyranges/reference/filter-ranges.md)
  :

  Subset a `Ranges` object

- [`arrange(`*`<Ranges>`*`)`](https://tidyomics.github.io/plyranges/reference/ranges-arrange.md)
  : Sort a Ranges object

- [`select(`*`<Ranges>`*`)`](https://tidyomics.github.io/plyranges/reference/ranges-select.md)
  : Select metadata columns of the Ranges object by name or position

- [`slice(`*`<Ranges>`*`)`](https://tidyomics.github.io/plyranges/reference/slice-ranges.md)
  [`slice(`*`<GroupedGenomicRanges>`*`)`](https://tidyomics.github.io/plyranges/reference/slice-ranges.md)
  [`slice(`*`<GroupedIntegerRanges>`*`)`](https://tidyomics.github.io/plyranges/reference/slice-ranges.md)
  : Choose rows by their position

- [`group_by(`*`<GenomicRanges>`*`)`](https://tidyomics.github.io/plyranges/reference/group_by-ranges.md)
  [`ungroup(`*`<GroupedGenomicRanges>`*`)`](https://tidyomics.github.io/plyranges/reference/group_by-ranges.md)
  [`groups(`*`<GroupedGenomicRanges>`*`)`](https://tidyomics.github.io/plyranges/reference/group_by-ranges.md)
  [`groups(`*`<GroupedIntegerRanges>`*`)`](https://tidyomics.github.io/plyranges/reference/group_by-ranges.md)
  : Group a Ranges by one or more variables

- [`pull(`*`<Ranges>`*`)`](https://tidyomics.github.io/plyranges/reference/pull-ranges.md)
  : Extract a single column from a Ranges object as a vector

- [`reduce_ranges()`](https://tidyomics.github.io/plyranges/reference/ranges-reduce.md)
  [`reduce_ranges_directed()`](https://tidyomics.github.io/plyranges/reference/ranges-reduce.md)
  : Reduce then aggregate a Ranges object

- [`disjoin_ranges()`](https://tidyomics.github.io/plyranges/reference/ranges-disjoin.md)
  [`disjoin_ranges_directed()`](https://tidyomics.github.io/plyranges/reference/ranges-disjoin.md)
  : Disjoin then aggregate a Ranges object

- [`compute_coverage()`](https://tidyomics.github.io/plyranges/reference/compute_coverage.md)
  : Compute coverage over a Ranges object

- [`expand_ranges()`](https://tidyomics.github.io/plyranges/reference/ranges-expand.md)
  : Expand list-columns in a Ranges object

## Overlaps

- [`join_overlap_intersect()`](https://tidyomics.github.io/plyranges/reference/overlap-joins.md)
  [`join_overlap_intersect_within()`](https://tidyomics.github.io/plyranges/reference/overlap-joins.md)
  [`join_overlap_intersect_directed()`](https://tidyomics.github.io/plyranges/reference/overlap-joins.md)
  [`join_overlap_intersect_within_directed()`](https://tidyomics.github.io/plyranges/reference/overlap-joins.md)
  [`join_overlap_inner()`](https://tidyomics.github.io/plyranges/reference/overlap-joins.md)
  [`join_overlap_inner_within()`](https://tidyomics.github.io/plyranges/reference/overlap-joins.md)
  [`join_overlap_inner_directed()`](https://tidyomics.github.io/plyranges/reference/overlap-joins.md)
  [`join_overlap_inner_within_directed()`](https://tidyomics.github.io/plyranges/reference/overlap-joins.md)
  [`join_overlap_left()`](https://tidyomics.github.io/plyranges/reference/overlap-joins.md)
  [`join_overlap_left_within()`](https://tidyomics.github.io/plyranges/reference/overlap-joins.md)
  [`join_overlap_left_directed()`](https://tidyomics.github.io/plyranges/reference/overlap-joins.md)
  [`join_overlap_left_within_directed()`](https://tidyomics.github.io/plyranges/reference/overlap-joins.md)
  : Join by overlapping Ranges
- [`join_overlap_self()`](https://tidyomics.github.io/plyranges/reference/ranges-overlaps-self.md)
  [`join_overlap_self_within()`](https://tidyomics.github.io/plyranges/reference/ranges-overlaps-self.md)
  [`join_overlap_self_directed()`](https://tidyomics.github.io/plyranges/reference/ranges-overlaps-self.md)
  [`join_overlap_self_within_directed()`](https://tidyomics.github.io/plyranges/reference/ranges-overlaps-self.md)
  : Find overlaps within a Ranges object
- [`filter_by_overlaps()`](https://tidyomics.github.io/plyranges/reference/ranges-filter-overlaps.md)
  [`filter_by_non_overlaps()`](https://tidyomics.github.io/plyranges/reference/ranges-filter-overlaps.md)
  [`filter_by_overlaps_directed()`](https://tidyomics.github.io/plyranges/reference/ranges-filter-overlaps.md)
  [`filter_by_non_overlaps_directed()`](https://tidyomics.github.io/plyranges/reference/ranges-filter-overlaps.md)
  : Filter by overlapping/non-overlapping ranges
- [`find_overlaps()`](https://tidyomics.github.io/plyranges/reference/ranges-overlaps.md)
  [`find_overlaps_within()`](https://tidyomics.github.io/plyranges/reference/ranges-overlaps.md)
  [`find_overlaps_directed()`](https://tidyomics.github.io/plyranges/reference/ranges-overlaps.md)
  [`find_overlaps_within_directed()`](https://tidyomics.github.io/plyranges/reference/ranges-overlaps.md)
  [`group_by_overlaps()`](https://tidyomics.github.io/plyranges/reference/ranges-overlaps.md)
  : Find overlap between two Ranges
- [`count_overlaps()`](https://tidyomics.github.io/plyranges/reference/ranges-count-overlaps.md)
  [`count_overlaps_within()`](https://tidyomics.github.io/plyranges/reference/ranges-count-overlaps.md)
  [`count_overlaps_directed()`](https://tidyomics.github.io/plyranges/reference/ranges-count-overlaps.md)
  [`count_overlaps_within_directed()`](https://tidyomics.github.io/plyranges/reference/ranges-count-overlaps.md)
  : Count the number of overlaps between two Ranges objects

## Joining on metadata columns

- [`join_mcols_inner()`](https://tidyomics.github.io/plyranges/reference/mcols-joins.md)
  [`join_mcols_left()`](https://tidyomics.github.io/plyranges/reference/mcols-joins.md)
  : Join data by metadata columns

## Nearest neighbours

- [`join_follow()`](https://tidyomics.github.io/plyranges/reference/ranges-follow.md)
  [`join_follow_left()`](https://tidyomics.github.io/plyranges/reference/ranges-follow.md)
  [`join_follow_upstream()`](https://tidyomics.github.io/plyranges/reference/ranges-follow.md)
  : Find following Ranges
- [`join_precede()`](https://tidyomics.github.io/plyranges/reference/ranges-precede.md)
  [`join_precede_right()`](https://tidyomics.github.io/plyranges/reference/ranges-precede.md)
  [`join_precede_downstream()`](https://tidyomics.github.io/plyranges/reference/ranges-precede.md)
  : Find preceding Ranges
- [`join_nearest()`](https://tidyomics.github.io/plyranges/reference/ranges-nearest.md)
  [`join_nearest_left()`](https://tidyomics.github.io/plyranges/reference/ranges-nearest.md)
  [`join_nearest_right()`](https://tidyomics.github.io/plyranges/reference/ranges-nearest.md)
  [`join_nearest_upstream()`](https://tidyomics.github.io/plyranges/reference/ranges-nearest.md)
  [`join_nearest_downstream()`](https://tidyomics.github.io/plyranges/reference/ranges-nearest.md)
  : Find nearest neighbours between two Ranges objects
- [`add_nearest_distance()`](https://tidyomics.github.io/plyranges/reference/add-nearest-distance.md)
  [`add_nearest_distance_left()`](https://tidyomics.github.io/plyranges/reference/add-nearest-distance.md)
  [`add_nearest_distance_right()`](https://tidyomics.github.io/plyranges/reference/add-nearest-distance.md)
  [`add_nearest_distance_upstream()`](https://tidyomics.github.io/plyranges/reference/add-nearest-distance.md)
  [`add_nearest_distance_downstream()`](https://tidyomics.github.io/plyranges/reference/add-nearest-distance.md)
  : Add distance to nearest neighbours between two Ranges objects

## Pairing up Ranges

- [`pair_overlaps()`](https://tidyomics.github.io/plyranges/reference/ranges-pairs.md)
  [`pair_nearest()`](https://tidyomics.github.io/plyranges/reference/ranges-pairs.md)
  [`pair_precede()`](https://tidyomics.github.io/plyranges/reference/ranges-pairs.md)
  [`pair_follow()`](https://tidyomics.github.io/plyranges/reference/ranges-pairs.md)
  : Pair together two ranges objects

## Set operations

- [`intersect_ranges()`](https://tidyomics.github.io/plyranges/reference/ranges-setops.md)
  [`intersect_ranges_directed()`](https://tidyomics.github.io/plyranges/reference/ranges-setops.md)
  [`union_ranges()`](https://tidyomics.github.io/plyranges/reference/ranges-setops.md)
  [`union_ranges_directed()`](https://tidyomics.github.io/plyranges/reference/ranges-setops.md)
  [`setdiff_ranges()`](https://tidyomics.github.io/plyranges/reference/ranges-setops.md)
  [`setdiff_ranges_directed()`](https://tidyomics.github.io/plyranges/reference/ranges-setops.md)
  [`complement_ranges()`](https://tidyomics.github.io/plyranges/reference/ranges-setops.md)
  [`complement_ranges_directed()`](https://tidyomics.github.io/plyranges/reference/ranges-setops.md)
  : Vector-wise Range set-operations
- [`` `%union%` ``](https://tidyomics.github.io/plyranges/reference/element-setops.md)
  [`` `%intersect%` ``](https://tidyomics.github.io/plyranges/reference/element-setops.md)
  [`` `%setdiff%` ``](https://tidyomics.github.io/plyranges/reference/element-setops.md)
  [`between()`](https://tidyomics.github.io/plyranges/reference/element-setops.md)
  [`span()`](https://tidyomics.github.io/plyranges/reference/element-setops.md)
  : Row-wise set operations on Ranges objects

## Annotation

- [`genome_info()`](https://tidyomics.github.io/plyranges/reference/ranges-info.md)
  [`set_genome_info()`](https://tidyomics.github.io/plyranges/reference/ranges-info.md)
  [`get_genome_info()`](https://tidyomics.github.io/plyranges/reference/ranges-info.md)
  : Construct annotation information

## Reading and Writing Genomics Data

- [`read_bam()`](https://tidyomics.github.io/plyranges/reference/io-bam-read.md)
  : Read a BAM file
- [`read_bed()`](https://tidyomics.github.io/plyranges/reference/io-bed-read.md)
  [`read_bed_graph()`](https://tidyomics.github.io/plyranges/reference/io-bed-read.md)
  [`read_narrowpeaks()`](https://tidyomics.github.io/plyranges/reference/io-bed-read.md)
  : Read a BED or BEDGraph file
- [`write_bed()`](https://tidyomics.github.io/plyranges/reference/io-bed-write.md)
  [`write_bed_graph()`](https://tidyomics.github.io/plyranges/reference/io-bed-write.md)
  [`write_narrowpeaks()`](https://tidyomics.github.io/plyranges/reference/io-bed-write.md)
  : Write a BED or BEDGraph file
- [`read_bigwig()`](https://tidyomics.github.io/plyranges/reference/io-bigwig-read.md)
  : Read a BigWig file
- [`write_bigwig()`](https://tidyomics.github.io/plyranges/reference/io-bigwig-write.md)
  : Write a BigWig file
- [`read_wig()`](https://tidyomics.github.io/plyranges/reference/io-wig-read.md)
  : Read a WIG file
- [`write_wig()`](https://tidyomics.github.io/plyranges/reference/io-wig-write.md)
  : Write a WIG file
- [`read_gff()`](https://tidyomics.github.io/plyranges/reference/io-gff-read.md)
  [`read_gff1()`](https://tidyomics.github.io/plyranges/reference/io-gff-read.md)
  [`read_gff2()`](https://tidyomics.github.io/plyranges/reference/io-gff-read.md)
  [`read_gff3()`](https://tidyomics.github.io/plyranges/reference/io-gff-read.md)
  : Read a GFF/GTF/GVT file
- [`write_gff()`](https://tidyomics.github.io/plyranges/reference/io-gff-write.md)
  [`write_gff1()`](https://tidyomics.github.io/plyranges/reference/io-gff-write.md)
  [`write_gff2()`](https://tidyomics.github.io/plyranges/reference/io-gff-write.md)
  [`write_gff3()`](https://tidyomics.github.io/plyranges/reference/io-gff-write.md)
  : Write a GFF(123) file
- [`FileOperator-class`](https://tidyomics.github.io/plyranges/reference/ranges-class.md)
  [`BamFileOperator-class`](https://tidyomics.github.io/plyranges/reference/ranges-class.md)
  : An abstract class to represent operations performed over a file
- [`DeferredGenomicRanges-class`](https://tidyomics.github.io/plyranges/reference/ranges-deferred.md)
  : DeferredGenomiRanges objects

## Utilities

- [`remove_names()`](https://tidyomics.github.io/plyranges/reference/ranges-names.md)
  [`names_to_column()`](https://tidyomics.github.io/plyranges/reference/ranges-names.md)
  [`id_to_column()`](https://tidyomics.github.io/plyranges/reference/ranges-names.md)
  : Tools for working with named Ranges
- [`set_width()`](https://tidyomics.github.io/plyranges/reference/ranges-setters.md)
  [`set_start()`](https://tidyomics.github.io/plyranges/reference/ranges-setters.md)
  [`set_end()`](https://tidyomics.github.io/plyranges/reference/ranges-setters.md)
  [`set_seqnames()`](https://tidyomics.github.io/plyranges/reference/ranges-setters.md)
  [`set_strand()`](https://tidyomics.github.io/plyranges/reference/ranges-setters.md)
  : Functional setters for Ranges objects
- [`as_ranges()`](https://tidyomics.github.io/plyranges/reference/as_ranges.md)
  : Coerce an Rle or RleList object to Ranges
- [`overscope_ranges()`](https://tidyomics.github.io/plyranges/reference/overscope_ranges.md)
  : Create an overscoped environment from a Ranges object
- [`` `%>%` ``](https://tidyomics.github.io/plyranges/reference/pipe.md)
  : Pipe operator
