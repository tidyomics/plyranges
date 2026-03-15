# Functional setters for Ranges objects

Functional setters for Ranges objects

## Usage

``` r
set_width(x, width)

set_start(x, start = 0L)

set_end(x, end = 0L)

set_seqnames(x, seqnames)

set_strand(x, strand)
```

## Arguments

- x:

  a Ranges object

- width:

  integer amount to modify width by

- start:

  integer amount to modify start by

- end:

  integer amount to modify end by

- seqnames:

  update seqnames column

- strand:

  update strand column

## Value

a Ranges object

## Details

These methods are used internally in `mutate()` to modify core columns
in Ranges objects.
