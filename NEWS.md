# plyranges 1.7.3

* patch left outer join when `x` or `y` are IRanges, flesh out overlaps documentation.

# plyranges 1.7.2

* Left outer join overlap operations now work if either `x` or `y` have
no metadata columns see [#70](https://github.com/sa-lee/plyranges/issues/70) 
* Left outer join overlap operations will also correctly
behave in situations when there are no non-overlapping ranges.
* Left outer join overlaps no longer modify seqinfo (see here)[https://support.bioconductor.org/p/125623/]

# plyranges 1.7.1

* Reformatting `NEWS.md` so no longer softlinks to inst/NEWS

# plyranges 1.5.13  

* plyranges release and devel have removed `unnest()` and replaced it
with `expand_ranges()` due to changes in the tidyr API. Please replace
all uses of this function with `expand_ranges()`

# plyranges 1.3.4 

* fixed bind_ranges so it preserves rownames

# plyranges 1.1.5 

* enable right generics to be called upon invoking plyranges functions
without loading plyranges

# plyranges 1.1.4 

* added tile/window methods
* fixed up documentation

# plyranges 1.1.3 

* doc updates

# plyranges 1.1.2 

* speed up of `group_by` methods
* refactor of BAM reading utilities

# plyranges 0.99.10 

* refactored `set_width` out so it's called internally by mutate
* along with `set_width` there are other internal `set_` methods
* add `_within_directed` variants for overlaps methods
* modified `overscope_ranges` to be an S3 method, should enable more refactoring in the future

# plyranges 0.99.9

https://bioconductor.org/packages/devel/bioc/html/plyranges.html

- package passed review and is now on Bioconductor devel branch!
- I've been pretty slack with updating the NEWS file but will be more
diligent in the future.