high to low priority

- [ ] gather tests from HelloRanges/IRanges/GRanges and port to plyranges
  - [ ] test-construction (making Ranges objects from df)
  - [ ] test-select
  - [ ] test-grouped classes
  - [ ] test-mutate
  - [ ] test-summarise
  - [ ] test-set-width
  - [ ] test-stretch
  - [ ] test-joins
  - [ ] test-nest
- [ ] include examples in docs
- [ ] vignette examples for all function suites
- [ ] `nranges` method for grouping
- [ ] find_overlaps should allow grouping
- [ ] overlaps directed methods
- [ ] sorting methods (i.e. `arrange`)
- [ ] update select method to use tidyselect::vars_select
- [ ] binary comparison operators
- [ ] add tidyselect to imports

capture `n()` call and port to length
look at aggregate methods 
test method for nest




DeferredGRanges class?

defer operations pushing down to data source -  avoid generics where possible

preserve semantics - if you mutate and add seqnames you should get back
a GRanges, dropping variables that break class should not be allowed 


