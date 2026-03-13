# Create an overscoped environment from a Ranges object

Create an overscoped environment from a Ranges object

## Usage

``` r
overscope_ranges(x, envir = parent.frame())
```

## Arguments

- x:

  a Ranges object

- envir:

  the environment to place the Ranges in (default =
  [`parent.frame()`](https://rdrr.io/r/base/sys.parent.html))

## Value

an environment

## Details

This is the backend for non-standard evaluation in `plyranges`.

## See also

[`rlang::new_data_mask()`](https://rlang.r-lib.org/reference/as_data_mask.html),
[`rlang::eval_tidy()`](https://rlang.r-lib.org/reference/eval_tidy.html)
