high to low priority

- [ ] implement sensible join operators - rename ojoin to *join_overlap to keep consitency
- [ ] update overlap join methods so that corresponding range information
is kept
- [ ] self_join_overlap function
- [ ] group_by_overlaps (group_by(ranges, query_hits)) and an `nranges` method for grouping
- [ ] find_overlaps should allow grouping
- [ ] nest/unest operators
- [ ] allow mutate to take grouped classes as input
- [ ] sorting methods (i.e. `arrange`)
- [ ] gather tests from HelloRanges/IRanges/GRanges and port to plyranges
- [ ] whole range set operators
- [ ] binary comparison operators
- [ ] vignette examples for all function suites




DeferredGRanges class?

defer operations pushing down to data source -  avoid generics where possible

preserve semantics - if you mutate and add seqnames you should get back
a GRanges, dropping variables that break class should not be allowed 


