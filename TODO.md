high to low priority

- [ ] group_by_overlaps (group_by(ranges, query_hits)) and an `nranges` method for grouping
- [ ] find_overlaps should allow grouping
- [ ] sorting methods (i.e. `arrange`)
- [ ] update select method to use tidyselect::vars_select
- [ ] gather tests from HelloRanges/IRanges/GRanges and port to plyranges
- [ ] binary comparison operators
- [ ] vignette examples for all function suites
- [ ] add tidyselect to imports



DeferredGRanges class?

defer operations pushing down to data source -  avoid generics where possible

preserve semantics - if you mutate and add seqnames you should get back
a GRanges, dropping variables that break class should not be allowed 


