# Group a Ranges by one or more variables

The function `group_by` takes a Ranges object and defines groups by one
or more variables. Operations are then performed on the Ranges by their
"group". `ungroup()` removes grouping.

## Usage

``` r
# S3 method for class 'GenomicRanges'
group_by(.data, ..., add = FALSE)

# S3 method for class 'GroupedGenomicRanges'
ungroup(x, ...)

# S3 method for class 'GroupedGenomicRanges'
groups(x)

# S3 method for class 'GroupedIntegerRanges'
groups(x)
```

## Arguments

- .data:

  a Ranges object.

- ...:

  Variable names to group by. These can be either metadata columns or
  the core variables of a Ranges.

- add:

  if `.data` is already a GroupedRanges object, when add = FALSE the
  (default), `group_by()` will override existing groups. If add = TRUE,
  additional groups will be added.

- x:

  a GroupedRanges object.

## Value

The `group_by()` function will return a GroupedRanges object. These have
the same appearance as a regular Ranges object but with an additional
groups slot.

## Details

`group_by()` creates a new object of class `GroupedGenomicRanges` if the
input is a `GRanges` object or an object of class `GroupedIntegerRanges`
if the input is a `IRanges` object. Both of these classes contain a slot
called `groups` corresponding to the names of grouping variables. They
also inherit from their parent classes, `Ranges` and `GenomicRanges`
respectively. `ungroup()` removes the grouping and will return either a
`GRanges` or `IRanges` object.

## Accessors

To return grouping variables on a grouped Ranges use either

- `groups(x)`:

  Returns a list of symbols

- `group_vars(x)`:

  Returns a character vector

## Examples

``` r
set.seed(100)
df <- data.frame(start = 1:10,
                 width = 5,
                 gc = runif(10),
                 cat = sample(letters[1:2], 10, replace = TRUE))
rng <- as_iranges(df)
rng_by_cat <- rng %>% group_by(cat)
# grouping does not change appearance or shape of Ranges
rng_by_cat
#> IRanges object with 10 ranges and 2 metadata columns:
#> Groups: cat [2]
#>            start       end     width |        gc         cat
#>        <integer> <integer> <integer> | <numeric> <character>
#>    [1]         1         5         5 | 0.3077661           b
#>    [2]         2         6         5 | 0.2576725           b
#>    [3]         3         7         5 | 0.5523224           b
#>    [4]         4         8         5 | 0.0563832           b
#>    [5]         5         9         5 | 0.4685493           a
#>    [6]         6        10         5 | 0.4837707           b
#>    [7]         7        11         5 | 0.8124026           b
#>    [8]         8        12         5 | 0.3703205           a
#>    [9]         9        13         5 | 0.5465586           a
#>   [10]        10        14         5 | 0.1702621           a
# a list of symbols
groups(rng_by_cat)
#> [[1]]
#> cat
#> 
# ungroup removes any grouping
ungroup(rng_by_cat)
#> IRanges object with 10 ranges and 2 metadata columns:
#>            start       end     width |        gc         cat
#>        <integer> <integer> <integer> | <numeric> <character>
#>    [1]         1         5         5 | 0.3077661           b
#>    [2]         2         6         5 | 0.2576725           b
#>    [3]         3         7         5 | 0.5523224           b
#>    [4]         4         8         5 | 0.0563832           b
#>    [5]         5         9         5 | 0.4685493           a
#>    [6]         6        10         5 | 0.4837707           b
#>    [7]         7        11         5 | 0.8124026           b
#>    [8]         8        12         5 | 0.3703205           a
#>    [9]         9        13         5 | 0.5465586           a
#>   [10]        10        14         5 | 0.1702621           a
# group_by works best with other verbs
grng <- as_granges(df,
                   seqnames = "chr1",
                   strand = sample(c("+", "-"), size = 10, replace = TRUE))

grng_by_strand <- grng %>% group_by(strand)
grng_by_strand
#> GRanges object with 10 ranges and 2 metadata columns:
#> Groups: strand [2]
#>        seqnames    ranges strand |        gc         cat
#>           <Rle> <IRanges>  <Rle> | <numeric> <character>
#>    [1]     chr1       1-5      + | 0.3077661           b
#>    [2]     chr1       2-6      - | 0.2576725           b
#>    [3]     chr1       3-7      - | 0.5523224           b
#>    [4]     chr1       4-8      + | 0.0563832           b
#>    [5]     chr1       5-9      - | 0.4685493           a
#>    [6]     chr1      6-10      + | 0.4837707           b
#>    [7]     chr1      7-11      + | 0.8124026           b
#>    [8]     chr1      8-12      - | 0.3703205           a
#>    [9]     chr1      9-13      - | 0.5465586           a
#>   [10]     chr1     10-14      + | 0.1702621           a
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
# grouping with other verbs
grng_by_strand %>% summarise(gc = mean(gc))
#> DataFrame with 2 rows and 2 columns
#>   strand        gc
#>    <Rle> <numeric>
#> 1      +  0.366117
#> 2      -  0.439085
grng_by_strand %>% filter(gc == min(gc))
#> GRanges object with 2 ranges and 2 metadata columns:
#> Groups: strand [2]
#>       seqnames    ranges strand |        gc         cat
#>          <Rle> <IRanges>  <Rle> | <numeric> <character>
#>   [1]     chr1       2-6      - | 0.2576725           b
#>   [2]     chr1       4-8      + | 0.0563832           b
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
grng_by_strand %>%
  ungroup() %>%
  summarise(gc = mean(gc))
#> DataFrame with 1 row and 1 column
#>          gc
#>   <numeric>
#> 1  0.402601

```
