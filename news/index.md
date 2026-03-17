# Changelog

## plyranges 1.31.5

- added `distance=TRUE` option for `join_overlap_*`
- added `join_mcols_left(x,y)` for GRanges `x` and tabular `y`

## plyranges 1.27.6

- added [`pull()`](https://dplyr.tidyverse.org/reference/pull.html)
  method for Ranges objects that extracts a single column as a vector
- moved `dplyr` to Depends and removing tidyverse reexports

## plyranges 1.9.3

- minor spelling and layout fixes to vignette,
- [@PeteHaitch](https://github.com/PeteHaitch) corrected table layout in
  vignette

## plyranges 1.9.2

- minor documentation fixes

## plyranges 1.9.1

- [Spencer Nystrom](https://github.com/snystrom) has made several
  significant contributions to the `join_nearest` family of functions:

  1.  `join_nearest_(x, y, ..., distance = TRUE)` family of functions
      now takes a new argument, `distance`, which allows the user to add
      a column for the distance of the nearest range `y` to that in `x`.
  2.  `add_nearest_distance_(x, y, ...)` family of functions, which will
      add a new metadata column to the `x` ranges object which contains
      the distance to its nearest neighbor in `y`. If there are no
      nearest neighbors, the new column will be given a missing value.

## plyranges 1.7.16

- refactoring of select internals, improved speed when casting a GRanges
  -\> DataFrame

## plyranges 1.7.15

- further fixes to reduce/disjoin internals

## plyranges 1.7.14

- fixes reduce/disjoin internals cleans up disjoin cases when an
  expansion occurs

## plyranges 1.7.13

- set tidyselect version to be v 1.0
- set coverage method for delegating ranges
- fix docs for bam reading

## plyranges 1.7.11

- move from
  [`tidyselect::vars_select()`](https://tidyselect.r-lib.org/reference/vars_select.html)
  to
  [`tidyselect::eval_select()`](https://tidyselect.r-lib.org/reference/eval_select.html)

## plyranges 1.7.7

- update handling of list columns,
  [`expand_ranges()`](https://tidyomics.github.io/plyranges/reference/ranges-expand.md)
  no longer takes cartesian product if lists are parallel.
  [`summarize()`](https://dplyr.tidyverse.org/reference/summarise.html)
  properly handles list column output without blowing out number of
  columns.

## plyranges 1.7.6

- adds method for
  [`dplyr::sample_n()`](https://dplyr.tidyverse.org/reference/sample_n.html)

## plyranges 1.7.5

- fixed issue [\#62](https://github.com/tidyomics/plyranges/issues/62)
  for Ranges construction, the
  [`as_granges()`](https://tidyomics.github.io/plyranges/reference/ranges-construct.md)
  and
  [`as_iranges()`](https://tidyomics.github.io/plyranges/reference/ranges-construct.md)
  functions now handle List columns correctly
- added in helper functions for dealing with names in Ranges. See
  `?ranges-names` for details.

## plyranges 1.7.4

- added [`slice()`](https://rdrr.io/pkg/IRanges/man/slice-methods.html)
  for Ranges, and GroupedRanges
- internals of grouping have been overhauled, but there shouldn’t be any
  user facing changes. It is now much faster to generate groupings.
- group information can be interrogated with
  [`dplyr::group_keys()`](https://dplyr.tidyverse.org/reference/group_data.html)
- a GRangesList can be obtained automatically from a
  GroupedGenomicRanges with
  [`dplyr::group_split()`](https://dplyr.tidyverse.org/reference/group_split.html)
- group indices can be generated with
  [`dplyr::group_indices()`](https://dplyr.tidyverse.org/reference/group_data.html)

## plyranges 1.7.3

- [`shift_downstream()`](https://tidyomics.github.io/plyranges/reference/ranges-shift.md)
  and
  [`shift_upstream()`](https://tidyomics.github.io/plyranges/reference/ranges-shift.md)
  now properly handle vector amounts of `shift`. Fixes issue
  [\#73](https://github.com/sa-lee/plyranges/issues/73)

## plyranges 1.7.2

- Left outer join overlap operations now work if either `x` or `y` have
  no metadata columns see
  [\#70](https://github.com/sa-lee/plyranges/issues/70)
- Left outer join overlap operations will also correctly behave in
  situations when there are no non-overlapping ranges.
- Left outer join overlaps no longer modify seqinfo (see
  here)\[<https://support.bioconductor.org/p/125623/>\]
- patch left outer join when `x` or `y` are IRanges, flesh out overlaps
  documentation.

## plyranges 1.7.1

- Reformatting `NEWS.md` so no longer softlinks to inst/NEWS

## plyranges 1.5.13

- plyranges release and devel have removed `unnest()` and replaced it
  with
  [`expand_ranges()`](https://tidyomics.github.io/plyranges/reference/ranges-expand.md)
  due to changes in the tidyr API. Please replace all uses of this
  function with
  [`expand_ranges()`](https://tidyomics.github.io/plyranges/reference/ranges-expand.md)

## plyranges 1.3.4

- fixed bind_ranges so it preserves rownames

## plyranges 1.1.5

- enable right generics to be called upon invoking plyranges functions
  without loading plyranges

## plyranges 1.1.4

- added tile/window methods
- fixed up documentation

## plyranges 1.1.3

- doc updates

## plyranges 1.1.2

- speed up of `group_by` methods
- refactor of BAM reading utilities

## plyranges 0.99.10

- refactored `set_width` out so it’s called internally by mutate
- along with `set_width` there are other internal `set_` methods
- add `_within_directed` variants for overlaps methods
- modified `overscope_ranges` to be an S3 method, should enable more
  refactoring in the future

## plyranges 0.99.9

<https://bioconductor.org/packages/devel/bioc/html/plyranges.html>

- package passed review and is now on Bioconductor devel branch!
- I’ve been pretty slack with updating the NEWS file but will be more
  diligent in the future.
