high to low priority

- [ ] gather tests from HelloRanges/IRanges/GRanges and port to plyranges
  - [ ] test-overlaps
  - [ ] test-nest
- [ ] include examples in docs for all exported functions
- [ ] fix errors and warnings in R CMD CHECK
- [ ] vignette examples for all function suites
- [ ] paper outline
- [ ] `nranges` method for grouping
- [ ] find_overlaps should allow grouping
- [ ] overlaps directed methods
- [ ] reduce/disjoin methods
- [ ] sorting methods (i.e. `arrange`)
- [ ] binary comparison operators
- [ ] add tidyselect to imports

capture `n()` call and port to length
look at aggregate methods 
test method for nest




DeferredGRanges class?

defer operations pushing down to data source -  avoid generics where possible

preserve semantics - if you mutate and add seqnames you should get back
a GRanges, dropping variables that break class should not be allowed 


